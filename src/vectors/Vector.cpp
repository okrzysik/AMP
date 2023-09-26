#include "AMP/vectors/Vector.h"
#include "AMP/IO/PIO.h"
#include "AMP/IO/RestartManager.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/VectorSelector.h"
#include "AMP/vectors/data/VectorDataNull.h"
#include "AMP/vectors/operations/VectorOperationsDefault.h"

#include <cfloat>
#include <cmath>
#include <cstring>
#include <typeinfo>


namespace AMP::LinearAlgebra {


/****************************************************************
 * Constructors/Destructor                                       *
 ****************************************************************/
Vector::Vector()
    : d_Variable( new Variable( "null" ) ),
      d_DOFManager( new AMP::Discretization::DOFManager( 0, AMP_MPI( AMP_COMM_SELF ) ) ),
      d_VectorData( new VectorDataNull<double>() ),
      d_VectorOps( new VectorOperationsDefault<double>() ),
      d_Views( new std::vector<std::any>() )
{
    AMPManager::incrementResource( "Vector" );
}
Vector::Vector( const std::string &name )
    : d_Variable( new Variable( name ) ),
      d_DOFManager( new AMP::Discretization::DOFManager( 0, AMP_MPI( AMP_COMM_SELF ) ) ),
      d_VectorData( new VectorDataNull<double>() ),
      d_VectorOps( new VectorOperationsDefault<double>() ),
      d_Views( new std::vector<std::any>() )
{
    AMPManager::incrementResource( "Vector" );
}
Vector::Vector( std::shared_ptr<VectorData> data,
                std::shared_ptr<VectorOperations> ops,
                std::shared_ptr<Variable> var,
                std::shared_ptr<AMP::Discretization::DOFManager> DOFManager )
    : d_Variable( var ),
      d_DOFManager( DOFManager ),
      d_VectorData( data ),
      d_VectorOps( ops ),
      d_Views( new std::vector<std::any>() )
{
    AMPManager::incrementResource( "Vector" );
    AMP_ASSERT( data && ops && var );
    if ( !d_DOFManager )
        d_DOFManager = std::make_shared<AMP::Discretization::DOFManager>(
            d_VectorData->getLocalSize(), d_VectorData->getComm() );
}
Vector::~Vector() { AMPManager::decrementResource( "Vector" ); }


/****************************************************************
 * Get the default name                                          *
 ****************************************************************/
std::string Vector::type() const { return "Vector<" + d_VectorData->VectorDataName() + ">"; }


/****************************************************************
 * setRandomValues                                               *
 ****************************************************************/
void Vector::setRandomValues() { d_VectorOps->setRandomValues( *getVectorData() ); }


/****************************************************************
 * Subset, View, and Select                                      *
 ****************************************************************/
Vector::shared_ptr Vector::selectInto( const VectorSelector &s )
{
    Vector::shared_ptr subvector;
    if ( s.isSelected( *this ) ) {
        // Subset the vector
        subvector = s.subset( shared_from_this() );
        if ( subvector ) {
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
    if ( s.isSelected( *this ) ) {
        // Subset the vector
        subvector = s.subset( shared_from_this() );
        if ( subvector ) {
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
        const std::string &name = dynamic_cast<const VS_ByVariableName *>( &s )->getName();
        if ( name == this->getVariable()->getName() )
            return shared_from_this();
    }
    Vector::shared_ptr retVal = this->selectInto( s );
    if ( retVal ) {
        if ( !std::dynamic_pointer_cast<MultiVector>( retVal ) )
            retVal = MultiVector::view( retVal, retVal->getComm() );
        auto var = std::make_shared<Variable>( variable_name );
        retVal->setVariable( var );
    }
    return retVal;
}
Vector::const_shared_ptr Vector::select( const VectorSelector &s,
                                         const std::string &variable_name ) const
{
    if ( dynamic_cast<const VS_ByVariableName *>( &s ) ) {
        const std::string &name = dynamic_cast<const VS_ByVariableName *>( &s )->getName();
        if ( name == this->getVariable()->getName() )
            return shared_from_this();
    }
    Vector::const_shared_ptr retVal = this->selectInto( s );
    if ( retVal ) {
        if ( !std::dynamic_pointer_cast<const MultiVector>( retVal ) )
            retVal = MultiVector::constView( retVal, retVal->getComm() );
        auto var = std::make_shared<Variable>( variable_name );
        std::const_pointer_cast<Vector>( retVal )->setVariable( var );
    }
    return retVal;
}
Vector::shared_ptr Vector::subsetVectorForVariable( const std::string &name )
{
    VS_ByVariableName selector( name );
    return selectInto( selector );
}
Vector::const_shared_ptr Vector::subsetVectorForVariable( const std::string &name ) const
{
    VS_ByVariableName selector( name );
    return selectInto( selector );
}

Vector::shared_ptr Vector::subsetVectorForVariable( std::shared_ptr<const Variable> var )
{
    if ( !var )
        return nullptr;
    auto selector = var->createVectorSelector();
    return selectInto( *selector );
}
Vector::const_shared_ptr
Vector::subsetVectorForVariable( std::shared_ptr<const Variable> var ) const
{
    if ( !var )
        return nullptr;
    auto selector = var->createVectorSelector();
    return selectInto( *selector );
}
Vector::shared_ptr Vector::subsetVectorForComponent( size_t index )
{
    return selectInto( VS_Components( index ) );
}
Vector::const_shared_ptr Vector::subsetVectorForComponent( size_t index ) const
{
    return selectInto( VS_Components( index ) );
}


/****************************************************************
 * clone, swap                                                   *
 ****************************************************************/
std::shared_ptr<Vector> Vector::clone() const { return clone( getVariable()->clone() ); }
std::shared_ptr<Vector> Vector::clone( const std::string &name ) const
{
    std::unique_ptr<Vector> retVal;
    if ( getVariable() ) {
        retVal = rawClone( getVariable()->clone( name ) );
    } else {
        retVal = rawClone( std::make_shared<Variable>( name ) );
    }
    return retVal;
}
Vector::shared_ptr Vector::clone( const std::shared_ptr<Variable> name ) const
{
    return rawClone( name );
}
std::unique_ptr<Vector> Vector::rawClone( const std::shared_ptr<Variable> name ) const
{
    auto vec          = std::make_unique<Vector>();
    vec->d_units      = d_units;
    vec->d_Variable   = name;
    vec->d_DOFManager = d_DOFManager;
    vec->d_VectorData = d_VectorData->cloneData();
    vec->d_VectorOps  = d_VectorOps->cloneOperations();
    if ( !vec->d_VectorData )
        AMP_ERROR( "Failed to clone data: " + d_VectorData->VectorDataName() );
    if ( !vec->d_VectorOps )
        AMP_ERROR( "Failed to clone ops: " + d_VectorOps->VectorOpName() );
    return vec;
}
void Vector::swapVectors( Vector &other )
{
    d_VectorData->swapData( *other.getVectorData() );
    std::swap( d_units, other.d_units );
}


/****************************************************************
 * Math API for Vector                                          *
 ****************************************************************/
void Vector::copy( const Vector &x ) { d_VectorOps->copy( *x.getVectorData(), *getVectorData() ); }
void Vector::zero() { d_VectorOps->zero( *getVectorData() ); }
void Vector::setToScalar( const Scalar &alpha )
{
    d_VectorOps->setToScalar( alpha, *getVectorData() );
}
void Vector::scale( const Scalar &alpha, const Vector &x )
{
    d_VectorOps->scale( alpha, *x.getVectorData(), *getVectorData() );
}
void Vector::scale( const Scalar &alpha ) { d_VectorOps->scale( alpha, *getVectorData() ); }
void Vector::add( const Vector &x, const Vector &y )
{
    d_VectorOps->add( *x.getVectorData(), *y.getVectorData(), *getVectorData() );
}
void Vector::subtract( const Vector &x, const Vector &y )
{
    d_VectorOps->subtract( *x.getVectorData(), *y.getVectorData(), *getVectorData() );
}
void Vector::multiply( const Vector &x, const Vector &y )
{
    d_VectorOps->multiply( *x.getVectorData(), *y.getVectorData(), *getVectorData() );
}
void Vector::divide( const Vector &x, const Vector &y )
{
    d_VectorOps->divide( *x.getVectorData(), *y.getVectorData(), *getVectorData() );
}
void Vector::reciprocal( const Vector &x )
{
    d_VectorOps->reciprocal( *x.getVectorData(), *getVectorData() );
}
void Vector::linearSum( const Scalar &alpha, const Vector &x, const Scalar &beta, const Vector &y )
{
    d_VectorOps->linearSum( alpha, *x.getVectorData(), beta, *y.getVectorData(), *getVectorData() );
}
void Vector::axpy( const Scalar &alpha, const Vector &x, const Vector &y )
{
    d_VectorOps->axpy( alpha, *x.getVectorData(), *y.getVectorData(), *getVectorData() );
}
void Vector::axpby( const Scalar &alpha, const Scalar &beta, const Vector &x )
{
    d_VectorOps->axpby( alpha, beta, *x.getVectorData(), *getVectorData() );
}
void Vector::abs( const Vector &x ) { d_VectorOps->abs( *x.getVectorData(), *getVectorData() ); }
void Vector::addScalar( const Vector &x, const Scalar &alpha_in )
{
    d_VectorOps->addScalar( *x.getVectorData(), alpha_in, *getVectorData() );
}
Scalar Vector::min() const { return d_VectorOps->min( *getVectorData() ); }
Scalar Vector::max() const { return d_VectorOps->max( *getVectorData() ); }
Scalar Vector::sum() const { return d_VectorOps->sum( *getVectorData() ); }
Scalar Vector::mean() const { return d_VectorOps->mean( *getVectorData() ); }
Scalar Vector::L1Norm() const { return d_VectorOps->L1Norm( *getVectorData() ); }
Scalar Vector::L2Norm() const { return d_VectorOps->L2Norm( *getVectorData() ); }
Scalar Vector::maxNorm() const { return d_VectorOps->maxNorm( *getVectorData() ); }
Scalar Vector::minQuotient( const Vector &x ) const
{
    return d_VectorOps->minQuotient( *x.getVectorData(), *getVectorData() );
}
Scalar Vector::wrmsNorm( const Vector &x, const Vector &y ) const
{
    return d_VectorOps->wrmsNorm( *x.getVectorData(), *y.getVectorData() );
}
Scalar Vector::wrmsNormMask( const Vector &x, const Vector &mask, const Vector &y ) const
{
    return d_VectorOps->wrmsNormMask(
        *x.getVectorData(), *mask.getVectorData(), *y.getVectorData() );
}
Scalar Vector::dot( const Vector &x ) const
{
    return d_VectorOps->dot( *getVectorData(), *x.getVectorData() );
}
bool Vector::equals( const Vector &a, const Scalar &tol ) const
{
    return d_VectorOps->equals( *a.getVectorData(), *getVectorData(), tol );
}


/****************************************************************
 * Get/set the name                                              *
 ****************************************************************/
std::string Vector::getName() const
{
    if ( d_Variable )
        return d_Variable->getName();
    return "";
}
void Vector::setName( const std::string &name )
{
    if ( d_Variable )
        d_Variable = d_Variable->clone( name );
    else
        d_Variable = std::make_shared<Variable>( name );
}

void Vector::rename( const std::string &src, const std::string &dst )
{
    if ( getName() == src )
        setName( dst );
    if ( dynamic_cast<MultiVector *>( this ) ) {
        for ( auto vec : *dynamic_cast<MultiVector *>( this ) )
            vec->rename( src, dst );
    }
}


/****************************************************************
 * Misc                                                          *
 ****************************************************************/
size_t Vector::getNumberOfComponents() const
{
    return d_VectorData ? d_VectorData->getNumberOfComponents() : 0;
}
uint64_t Vector::getID() const { return getComm().bcast( reinterpret_cast<uint64_t>( this ), 0 ); }
void Vector::setUnits( AMP::Units units ) { d_units = units; }
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


} // namespace AMP::LinearAlgebra


/********************************************************
 *  Restart operations                                   *
 ********************************************************/
template<>
AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::Vector>::DataStoreType(
    const std::string &name,
    std::shared_ptr<const AMP::LinearAlgebra::Vector> vec,
    RestartManager *manager )
    : d_data( vec )
{
    d_name = name;
    d_hash = vec->getID();
    manager->registerData( vec->getVariable() );
    manager->registerData( vec->getDOFManager() );
    manager->registerData( vec->getVectorData() );
    manager->registerData( vec->getVectorOperations() );
    auto multivec = std::dynamic_pointer_cast<const AMP::LinearAlgebra::MultiVector>( vec );
    if ( multivec ) {
        for ( auto vec2 : *multivec )
            manager->registerData( vec2 );
    }
}
template<>
void AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::Vector>::write(
    hid_t fid, const std::string &name ) const
{
    hid_t gid = createGroup( fid, name );
    writeHDF5( gid, "units", d_data->getUnits() );
    writeHDF5( gid, "var", d_data->getVariable()->getID() );
    writeHDF5( gid, "dofs", d_data->getDOFManager()->getID() );
    writeHDF5( gid, "data", d_data->getVectorData()->getID() );
    writeHDF5( gid, "ops", d_data->getVectorOperations()->getID() );
    auto multivec = std::dynamic_pointer_cast<const AMP::LinearAlgebra::MultiVector>( d_data );
    if ( multivec ) {
        std::vector<uint64_t> vecs;
        for ( auto vec2 : *multivec )
            vecs.push_back( vec2->getID() );
        writeHDF5( gid, "vecs", vecs );
    }
    closeGroup( gid );
}
template<>
std::shared_ptr<AMP::LinearAlgebra::Vector>
AMP::IO::RestartManager::getData<AMP::LinearAlgebra::Vector>( uint64_t hash )
{
    hid_t gid = openGroup( d_fid, hash2String( hash ) );
    AMP_ERROR( "Not finished" );
    closeGroup( gid );
}
