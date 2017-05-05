#include "vectors/Vector.h"
#include "utils/AMP_MPI.h"
#include "utils/Counter.h"
#include "utils/Utilities.h"
#include "vectors/DataChangeFirer.h"
#include "vectors/MultiVector.h"
#include "vectors/VectorSelector.h"

#include <float.h>
#include <math.h>
#include <string.h>
#include <typeinfo>


namespace AMP {
namespace LinearAlgebra {


RNG::shared_ptr Vector::d_DefaultRNG;
#define DESCRIPTOR_ID_ARRAY_SCRATCH_SPACE ( 10 )


/****************************************************************
* Constructors                                                  *
****************************************************************/
Vector::Vector():
    VectorData(),
    VectorOperations()
{
    d_VectorData = dynamic_cast<VectorData*>(this);
    d_Ghosts    = AMP::shared_ptr<std::vector<double>>( new std::vector<double> );
    d_AddBuffer = AMP::shared_ptr<std::vector<double>>( new std::vector<double> );
    d_UpdateState.reset( new UpdateState );
    *d_UpdateState = UpdateState::UNCHANGED;
    d_Views        = AMP::shared_ptr<std::vector<AMP::weak_ptr<Vector>>>(
        new std::vector<AMP::weak_ptr<Vector>>() );
}
Vector::Vector( VectorParameters::shared_ptr parameters )
{
    d_VectorData = dynamic_cast<VectorData*>(this);
    // Set default output stream
    d_output_stream = &AMP::plog;
    // Copy the relavent parameters
    AMP_INSIST( parameters->d_CommList.get() != nullptr,
                "d_CommList must be set in VectorParameters" );
    AMP_INSIST( parameters->d_DOFManager.get() != nullptr,
                "d_DOFManager must be set in VectorParameters" );
    setCommunicationList( parameters->d_CommList );
    d_UpdateState.reset( new UpdateState );
    *d_UpdateState = UpdateState::UNCHANGED;
    d_DOFManager   = parameters->d_DOFManager;
    d_Views        = AMP::shared_ptr<std::vector<AMP::weak_ptr<Vector>>>(
        new std::vector<AMP::weak_ptr<Vector>>() );
}


/****************************************************************
* De-Constructors                                               *
****************************************************************/
Vector::~Vector() {}


/****************************************************************
* Subset, View, and Select                                      *
****************************************************************/
Vector::shared_ptr Vector::selectInto( const VectorSelector &s )
{
    Vector::shared_ptr subvector;
    if ( s.isSelected( shared_from_this() ) ) {
        // Subset the vector
        subvector = s.subset( shared_from_this() );
        if ( subvector != nullptr ) {
            // Check the global size of the new vector to make sure it is <= the current size
            size_t N1 = this->getGlobalSize();
            size_t N2 = subvector->getGlobalSize();
            AMP_ASSERT( N2 <= N1 );
        }
    }
    return subvector;
}
Vector::const_shared_ptr Vector::selectInto( const VectorSelector &s ) const
{
    Vector::const_shared_ptr subvector;
    if ( s.isSelected( shared_from_this() ) ) {
        // Subset the vector
        subvector = s.subset( shared_from_this() );
        if ( subvector != nullptr ) {
            // Check the global size of the new vector to make sure it is <= the current size
            size_t N1 = this->getGlobalSize();
            size_t N2 = subvector->getGlobalSize();
            AMP_ASSERT( N2 <= N1 );
        }
    }
    return subvector;
}
Vector::shared_ptr Vector::select( const VectorSelector &s, const std::string &variable_name )
{
    if ( dynamic_cast<const VS_ByVariableName *>( &s ) ) {
        std::string name = dynamic_cast<const VS_ByVariableName *>( &s )->getName();
        if ( name == this->getVariable()->getName() )
            return shared_from_this();
    }
    Vector::shared_ptr retVal = this->selectInto( s );
    if ( retVal != nullptr ) {
        if ( AMP::dynamic_pointer_cast<MultiVector>( retVal ) == nullptr )
            retVal = MultiVector::view( retVal, retVal->getComm() );
        Variable::shared_ptr var( new Variable( variable_name ) );
        retVal->setVariable( var );
    }
    return retVal;
}
Vector::const_shared_ptr Vector::constSelect( const VectorSelector &s,
                                              const std::string &variable_name ) const
{
    if ( dynamic_cast<const VS_ByVariableName *>( &s ) ) {
        std::string name = dynamic_cast<const VS_ByVariableName *>( &s )->getName();
        if ( name == this->getVariable()->getName() )
            return shared_from_this();
    }
    Vector::const_shared_ptr retVal = this->selectInto( s );
    if ( retVal != nullptr ) {
        if ( AMP::dynamic_pointer_cast<const MultiVector>( retVal ) == nullptr )
            retVal = MultiVector::view( retVal, retVal->getComm() );
        Variable::shared_ptr var( new Variable( variable_name ) );
        AMP::const_pointer_cast<Vector>( retVal )->setVariable( var );
    }
    return retVal;
}
void Vector::registerView( Vector::shared_ptr v ) const
{
    for ( size_t i = 0; i != d_Views->size(); i++ )
        if ( ( *d_Views )[i].lock() == v )
            return;
    ( *d_Views ).push_back( v );
}
Vector::shared_ptr Vector::subsetVectorForVariable( const Variable::shared_ptr &name )
{
    Vector::shared_ptr retVal;
    if ( d_pVariable ) { // If there is a variable...
        if ( *d_pVariable == *name )
            retVal = shared_from_this();
    }
    return retVal;
}
Vector::const_shared_ptr
Vector::constSubsetVectorForVariable( const Variable::shared_ptr &name ) const
{
    Vector::const_shared_ptr retVal;
    if ( d_pVariable ) { // If there is a variable...
        if ( *d_pVariable == *name )
            retVal = shared_from_this();
    }
    return retVal;
}
Vector::shared_ptr Vector::cloneVector( const std::string &name ) const
{
    Vector::shared_ptr retVal;
    if ( getVariable() ) {
        retVal = cloneVector( getVariable()->cloneVariable( name ) );
    } else {
        retVal = cloneVector( Variable::shared_ptr( new Variable( name ) ) );
    }
    return retVal;
}


/****************************************************************
* Misc                                                          *
****************************************************************/
void Vector::copyVector( Vector::const_shared_ptr rhs )
{
    if ( rhs->getLocalSize() != getLocalSize() )
        AMP_ERROR( "Destination vector and source vector not the same size" );
    VectorDataIterator cur1      = begin();
    VectorDataIterator end1      = end();
    ConstVectorDataIterator cur2 = rhs->begin();
    while ( cur1 != end1 ) {
        *cur1 = *cur2;
        ++cur1;
        ++cur2;
    }
    if ( isA<DataChangeFirer>() )
        castTo<DataChangeFirer>().fireDataChange();
    copyGhostValues( rhs );
    // Copy the consistency state from the rhs
    *d_UpdateState = *( rhs->getUpdateStatusPtr() );
}
void Vector::setCommunicationList( CommunicationList::shared_ptr comm )
{
    AMP_ASSERT( comm );
    d_CommList = comm;
    if ( comm ) {
        addCommunicationListToParameters( comm );
        d_Ghosts = AMP::shared_ptr<std::vector<double>>(
            new std::vector<double>( d_CommList->getVectorReceiveBufferSize() ) );
        d_AddBuffer = AMP::shared_ptr<std::vector<double>>(
            new std::vector<double>( d_CommList->getVectorReceiveBufferSize() ) );
    }
}
bool Vector::equals( Vector const &rhs, double tol ) const
{
    int RetVal = 0;
    if ( ( getGlobalSize() == rhs.getGlobalSize() ) && ( getLocalSize() == rhs.getLocalSize() ) ) {
        ConstVectorDataIterator cur1 = begin();
        ConstVectorDataIterator cur2 = rhs.begin();
        ConstVectorDataIterator last = end();
        bool failed                  = false;
        while ( cur1 != last ) {
            double v1 = *cur1;
            double v2 = *cur2;
            if ( fabs( v1 - v2 ) > tol ) {
                failed = true;
                break;
            }
            ++cur1;
            ++cur2;
        }
        if ( !failed ) {
            RetVal = 1;
        }
    }

    int ans;
    if ( d_CommList ) {
        ans = d_CommList->getComm().minReduce( RetVal );
    } else {
        ans = RetVal;
    }

    return ans == 1 ? true : false;
}
void Vector::copyGhostValues( const AMP::shared_ptr<const Vector> &rhs )
{
    if ( getGhostSize() == 0 ) {
        // No ghosts to fill, we don't need to do anything
    } else if ( getGhostSize() == rhs->getGhostSize() ) {
        // The ghosts in the src vector match the current vector
        // Copy the ghosts from the rhs
        std::vector<size_t> ghostIDs = getCommunicationList()->getGhostIDList();
        std::vector<double> values( ghostIDs.size() );
        rhs->getGhostValuesByGlobalID( ghostIDs.size(), &ghostIDs[0], &values[0] );
        this->setGhostValuesByGlobalID( ghostIDs.size(), &ghostIDs[0], &values[0] );
        // Copy the consistency state from the rhs
        *d_UpdateState = *( rhs->getUpdateStatusPtr() );
    } else {
        // We can't copy the ghosts from the rhs
        // Use makeConsistent to fill the ghosts
        // Note: this will incure global communication
        makeConsistent( ScatterType::CONSISTENT_SET );
    }
}


std::ostream &operator<<( std::ostream &out, const Vector &v )
{
    out << "Vector type: " << v.type() << "\n";
    if ( v.getVariable() ) {
        out << "Variable name: " << v.getVariable()->getName() << "\n";
    }
    if ( v.getCommunicationList() ) {
        int rank = v.getComm().getRank();
        out << "Processor: " << rank << "\n";
    }
    out << "\n"
        << "Number of owned elements: " << v.getLocalSize() << "\n"
        << "Number of ghosted elements: " << v.getGhostSize() << "\n";
    if ( v.numberOfDataBlocks() > 1 )
        out << "Number of sub-vectors: " << v.numberOfDataBlocks() << "\n";
    out << "\n"
        << "Data block pointers: \n";

    for ( size_t i = 0; i != v.numberOfDataBlocks(); i++ )
        out << "                     " << v.getRawDataBlock<double>( i ) << "  ( "
            << v.sizeOfDataBlock( i ) << " elements )\n";

    out << "\n"
        << "Local Data:\n";
    v.dumpOwnedData( out );

    out << "\n"
        << "Ghosted Data:\n";
    v.dumpGhostedData( out );

    return out;
}


} // LinearAlgebra namespace
} // AMP namespace

