#include "AMP/vectors/Vector.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/VectorSelector.h"
#include "AMP/vectors/data/VectorDataNull.h"
#include "AMP/vectors/operations/VectorOperationsDefault.h"

#include <cfloat>
#include <cmath>
#include <cstring>
#include <typeinfo>


namespace AMP {
namespace LinearAlgebra {


RNG::shared_ptr Vector::d_DefaultRNG;


/****************************************************************
 * Constructors/Destructor                                       *
 ****************************************************************/
Vector::Vector()
    : d_pVariable( new Variable( "null" ) ),
      d_DOFManager( new AMP::Discretization::DOFManager( 0, AMP_MPI( AMP_COMM_SELF ) ) ),
      d_VectorData( new VectorDataNull<double>() ),
      d_VectorOps( new VectorOperationsDefault<double>() ),
      d_Views( new std::vector<std::weak_ptr<Vector>>() ),
      d_output_stream( &AMP::plog )
{
}
Vector::Vector( const std::string &name )
    : d_pVariable( new Variable( name ) ),
      d_DOFManager( new AMP::Discretization::DOFManager( 0, AMP_MPI( AMP_COMM_SELF ) ) ),
      d_VectorData( new VectorDataNull<double>() ),
      d_VectorOps( new VectorOperationsDefault<double>() ),
      d_Views( new std::vector<std::weak_ptr<Vector>>() ),
      d_output_stream( &AMP::plog )
{
}
Vector::Vector( std::shared_ptr<VectorData> data,
                std::shared_ptr<VectorOperations> ops,
                Variable::shared_ptr var,
                AMP::Discretization::DOFManager::shared_ptr DOFManager )
    : d_pVariable( var ),
      d_DOFManager( DOFManager ),
      d_VectorData( data ),
      d_VectorOps( ops ),
      d_Views( new std::vector<std::weak_ptr<Vector>>() ),
      d_output_stream( &AMP::plog )
{
    AMP_ASSERT( data && ops && var );
    if ( !d_DOFManager )
        d_DOFManager = std::make_shared<AMP::Discretization::DOFManager>(
            d_VectorData->getLocalSize(), d_VectorData->getComm() );
}
Vector::~Vector() {}


/****************************************************************
 * Get the default name                                          *
 ****************************************************************/
std::string Vector::type() const { return "Vector<" + d_VectorData->VectorDataName() + ">"; }


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
    if ( !retVal )
        printf( "Unable to subset for %s in %s\n",
                name->getName().data(),
                getVariable()->getName().data() );
    return retVal;
}


/****************************************************************
 * clone, swap                                                   *
 ****************************************************************/
Vector::shared_ptr Vector::cloneVector() const { return cloneVector( getVariable() ); }
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
Vector::shared_ptr Vector::cloneVector( const Variable::shared_ptr name ) const
{
    auto vec             = std::make_shared<Vector>();
    vec->d_pVariable     = name;
    vec->d_DOFManager    = d_DOFManager;
    vec->d_VectorData    = d_VectorData->cloneData();
    vec->d_VectorOps     = d_VectorOps->cloneOperations();
    vec->d_output_stream = d_output_stream;
    return vec;
}
void Vector::swapVectors( Vector &other ) { d_VectorData->swapData( *other.getVectorData() ); }


/****************************************************************
 * Misc                                                          *
 ****************************************************************/
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
