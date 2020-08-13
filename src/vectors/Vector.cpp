#include "AMP/vectors/Vector.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/DataChangeFirer.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/VectorSelector.h"
#include "AMP/vectors/operations/VectorOperationsDefault.h"

#include <cfloat>
#include <cmath>
#include <cstring>
#include <typeinfo>


namespace AMP {
namespace LinearAlgebra {


RNG::shared_ptr Vector::d_DefaultRNG;
#define DESCRIPTOR_ID_ARRAY_SCRATCH_SPACE ( 10 )


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
Vector::Vector() : VectorData()
{
    d_VectorData = dynamic_cast<VectorData *>( this );
    d_Ghosts     = std::make_shared<std::vector<double>>();
    d_AddBuffer  = std::make_shared<std::vector<double>>();
    d_UpdateState.reset( new UpdateState );
    *d_UpdateState = UpdateState::UNCHANGED;
    d_Views        = std::make_shared<std::vector<std::weak_ptr<Vector>>>();
    // Set default output stream
    d_output_stream = &AMP::plog;
}
Vector::Vector( VectorParameters::shared_ptr parameters )
{
    d_VectorData = dynamic_cast<VectorData *>( this );
    // Set default output stream
    d_output_stream = &AMP::plog;
    // Copy the relavent parameters
    AMP_INSIST( parameters->d_CommList, "d_CommList must be set in VectorParameters" );
    AMP_INSIST( parameters->d_DOFManager, "d_DOFManager must be set in VectorParameters" );
    setCommunicationList( parameters->d_CommList );
    d_UpdateState.reset( new UpdateState );
    *d_UpdateState = UpdateState::UNCHANGED;
    d_DOFManager   = parameters->d_DOFManager;
    d_Views        = std::make_shared<std::vector<std::weak_ptr<Vector>>>();
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
        if ( std::dynamic_pointer_cast<MultiVector>( retVal ) == nullptr )
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
        if ( std::dynamic_pointer_cast<const MultiVector>( retVal ) == nullptr )
            retVal = MultiVector::constView( retVal, retVal->getComm() );
        Variable::shared_ptr var( new Variable( variable_name ) );
        std::const_pointer_cast<Vector>( retVal )->setVariable( var );
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
Vector::shared_ptr Vector::subsetVectorForVariable( Variable::const_shared_ptr name )
{
    Vector::shared_ptr retVal;
    if ( d_pVariable ) { // If there is a variable...
        if ( *d_pVariable == *name )
            retVal = shared_from_this();
    }
    return retVal;
}
Vector::const_shared_ptr
Vector::constSubsetVectorForVariable( Variable::const_shared_ptr name ) const
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
        retVal = cloneVector( std::make_shared<Variable>( name ) );
    }
    return retVal;
}


/****************************************************************
 * Misc                                                          *
 ****************************************************************/
void Vector::setCommunicationList( CommunicationList::shared_ptr comm )
{
    AMP_ASSERT( comm );
    d_CommList = comm;
    if ( comm ) {
        addCommunicationListToParameters( comm );
        d_Ghosts =
            std::make_shared<std::vector<double>>( d_CommList->getVectorReceiveBufferSize() );
        d_AddBuffer =
            std::make_shared<std::vector<double>>( d_CommList->getVectorReceiveBufferSize() );
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


} // namespace LinearAlgebra
} // namespace AMP
