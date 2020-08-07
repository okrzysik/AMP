#include "AMP/vectors/trilinos/thyra/NativeThyraVector.h"
#include "AMP/vectors/trilinos/thyra/ThyraVectorWrapper.h"


// Trilinos includes
DISABLE_WARNINGS
//#include "Thyra_SpmdVectorBase_def.hpp"
#include "Thyra_DefaultSpmdVector_def.hpp"
#include "Thyra_VectorStdOps_def.hpp"
#include "Trilinos_version.h"
ENABLE_WARNINGS


namespace AMP {
namespace LinearAlgebra {


/************************************************************************
 * Constructors                                                          *
 ************************************************************************/
NativeThyraVector::NativeThyraVector( VectorParameters::shared_ptr in_params )
    : Vector(), ThyraVector()
{
    auto params = std::dynamic_pointer_cast<NativeThyraVectorParameters>( in_params );
    AMP_ASSERT( params != nullptr );
    AMP_ASSERT( !params->d_comm.isNull() );
    AMP_ASSERT( params->d_InVec.get() != nullptr );
    Thyra::Ordinal dim = params->d_InVec->space()->dim();
    AMP_ASSERT( params->d_comm.sumReduce( params->d_local ) == static_cast<size_t>( dim ) );
    auto communicationListParams         = std::make_shared<CommunicationListParameters>();
    communicationListParams->d_comm      = params->d_comm;
    communicationListParams->d_localsize = params->d_local;
    d_CommList   = std::make_shared<CommunicationList>( communicationListParams );
    d_DOFManager = std::make_shared<Discretization::DOFManager>( params->d_local, params->d_comm );
    d_local      = params->d_local;
    d_thyraVec   = params->d_InVec;
    d_pVariable  = params->d_var;
}


/************************************************************************
 * Destructor                                                            *
 ************************************************************************/
NativeThyraVector::~NativeThyraVector() = default;


/************************************************************************
 * Vector functions                                                      *
 ************************************************************************/
Vector::shared_ptr NativeThyraVector::cloneVector( const Variable::shared_ptr var ) const
{
    std::shared_ptr<NativeThyraVectorParameters> params( new NativeThyraVectorParameters() );
    params->d_InVec = d_thyraVec->clone_v();
    params->d_local = d_local;
    params->d_comm  = getComm();
    params->d_var   = var;
    return std::make_shared<NativeThyraVector>( params );
}


void NativeThyraVector::putRawData( const double *in )
{
    size_t i = 0;
    for ( size_t b = 0; b < numberOfDataBlocks(); b++ ) {
        auto *data = reinterpret_cast<double *>( getRawDataBlockAsVoid( b ) );
        for ( size_t j = 0; j < sizeOfDataBlock( b ); j++, i++ )
            data[j] = in[i];
    }
}


void NativeThyraVector::copyOutRawData( double *out ) const
{
    size_t i = 0;
    for ( size_t b = 0; b < numberOfDataBlocks(); b++ ) {
        const auto *data = reinterpret_cast<const double *>( getRawDataBlockAsVoid( b ) );
        for ( size_t j = 0; j < sizeOfDataBlock( b ); j++, i++ )
            out[i] = data[j];
    }
}


void *NativeThyraVector::getRawDataBlockAsVoid( size_t i )
{
    Thyra::VectorBase<double> *ptr = d_thyraVec.get();
    auto *spmdVector               = dynamic_cast<Thyra::DefaultSpmdVector<double> *>( ptr );
    if ( spmdVector != nullptr ) {
        if ( i != 0 )
            AMP_ERROR( "Invalid block" );
        return spmdVector->getPtr();
    }
    auto *wrapperVector = dynamic_cast<ThyraVectorWrapper *>( ptr );
    if ( wrapperVector != nullptr ) {
        AMP_INSIST( wrapperVector->numVecs() == 1,
                    "Not ready for dealing with multiple copies of the vector yet" );
        return wrapperVector->getVec( 0 )->getRawDataBlock<double>( i );
    }
    AMP_ERROR( "not finished" );
    return nullptr;
}


const void *NativeThyraVector::getRawDataBlockAsVoid( size_t i ) const
{
    const Thyra::VectorBase<double> *ptr = d_thyraVec.get();
    const auto *spmdVector = dynamic_cast<const Thyra::DefaultSpmdVector<double> *>( ptr );
    if ( spmdVector != nullptr ) {
        if ( i != 0 )
            AMP_ERROR( "Invalid block" );
        return spmdVector->getPtr();
    }
    const auto *wrapperVector = dynamic_cast<const ThyraVectorWrapper *>( ptr );
    if ( wrapperVector != nullptr ) {
        AMP_INSIST( wrapperVector->numVecs() == 1,
                    "Not ready for dealing with multiple copies of the vector yet" );
        return wrapperVector->getVec( 0 )->getRawDataBlock<double>( i );
    }
    return nullptr;
}


size_t NativeThyraVector::numberOfDataBlocks() const
{
    const Thyra::VectorBase<double> *ptr = d_thyraVec.get();
    if ( dynamic_cast<const Thyra::DefaultSpmdVector<double> *>( ptr ) != nullptr )
        return 1;
    const auto *wrapperVector = dynamic_cast<const ThyraVectorWrapper *>( ptr );
    if ( wrapperVector != nullptr )
        return wrapperVector->getVec( 0 )->numberOfDataBlocks();
    AMP_ERROR( "not finished" );
    return 1;
}


size_t NativeThyraVector::sizeOfDataBlock( size_t i ) const
{
    const Thyra::VectorBase<double> *ptr = d_thyraVec.get();
    const auto *spmdVector = dynamic_cast<const Thyra::DefaultSpmdVector<double> *>( ptr );
    if ( spmdVector != nullptr )
        return d_local;
    const auto *wrapperVector = dynamic_cast<const ThyraVectorWrapper *>( ptr );
    if ( wrapperVector != nullptr ) {
        AMP_INSIST( wrapperVector->numVecs() == 1,
                    "Not ready for dealing with multiple copies of the vector yet" );
        return wrapperVector->getVec( 0 )->sizeOfDataBlock( i );
    }
    AMP_ERROR( "not finished" );
    return d_local;
}

void NativeThyraVector::getValuesByLocalID( int numVals, size_t *ndx, double *vals ) const
{
    Vector::getValuesByLocalID( numVals, ndx, vals );
}

void NativeThyraVector::aliasVector( Vector & ) { AMP_ERROR( "not implemented" ); }


void NativeThyraVector::swapVectors( Vector & ) { AMP_ERROR( "not implemented" ); }


Teuchos::RCP<const Thyra::VectorBase<double>>
NativeThyraVector::getThyraVec( const VectorOperations &v )
{
    auto vec = dynamic_cast<const Vector *>( &v );
    AMP_ASSERT( vec != nullptr );
    auto vec2 = std::dynamic_pointer_cast<const ThyraVector>(
        ThyraVector::constView( vec->shared_from_this() ) );
    AMP_ASSERT( vec2 != nullptr );
    return vec2->getVec();
}


Teuchos::RCP<const Thyra::VectorBase<double>>
NativeThyraVector::getThyraVec( const Vector::const_shared_ptr &vec )
{
    auto vec2 = std::dynamic_pointer_cast<const ThyraVector>( ThyraVector::constView( vec ) );
    AMP_ASSERT( vec2 != nullptr );
    return vec2->getVec();
}

Teuchos::RCP<const Thyra::VectorBase<double>>
NativeThyraVector::getThyraVec( const VectorData &v )
{
    auto vec = dynamic_cast<const Vector *>( &v );
    AMP_ASSERT( vec != nullptr );
    auto vec2 = std::dynamic_pointer_cast<const ThyraVector>(
        ThyraVector::constView( vec->shared_from_this() ) );
    AMP_ASSERT( vec2 != nullptr );
    return vec2->getVec();
}

Teuchos::RCP<Thyra::VectorBase<double>>
NativeThyraVector::getThyraVec( VectorData &v )
{
    auto vec = dynamic_cast<Vector *>( &v );
    AMP_ASSERT( vec != nullptr );
    auto vec2 = std::dynamic_pointer_cast<ThyraVector>(
        ThyraVector::view( vec->shared_from_this() ) );
    AMP_ASSERT( vec2 != nullptr );
    return vec2->getVec();
}
#if 0
void NativeThyraVector::copy( const VectorOperations &src_vec )
{
  copy( *(src_vec.getVectorData()), *getVectorData() );
}

void NativeThyraVector::setToScalar( double alpha )
{
  setToScalar(alpha, *getVectorData() );
}

void NativeThyraVector::setRandomValues( void )
{
   setRandomValues( *getVectorData() );
}

void NativeThyraVector::scale( double alpha, const VectorOperations &x )
{
  scale(alpha, *(x.getVectorData()), *getVectorData());
}


void NativeThyraVector::scale( double alpha )
{
  scale(alpha, *getVectorData());
}

void NativeThyraVector::add( const VectorOperations &x, const VectorOperations &y )
{
  add( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}


void NativeThyraVector::subtract( const VectorOperations &x, const VectorOperations &y )
{
  subtract( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}


void NativeThyraVector::multiply( const VectorOperations &x, const VectorOperations &y )
{
  multiply( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}


void NativeThyraVector::divide( const VectorOperations &x, const VectorOperations &y )
{
  divide( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}


void NativeThyraVector::reciprocal( const VectorOperations &x )
{
  reciprocal( *(x.getVectorData()), *getVectorData() );
}


void NativeThyraVector::linearSum( double alpha,
                                   const VectorOperations &x,
                                   double beta,
                                   const VectorOperations &y )
{
  linearSum( alpha,
	     *(x.getVectorData()),
	     beta,
	     *(y.getVectorData()),
	     *getVectorData() );
}


void NativeThyraVector::axpy( double alpha, const VectorOperations &x, const VectorOperations &y )
{
  axpy( alpha,
	*(x.getVectorData()),
	*(y.getVectorData()),
	*getVectorData() );
}


void NativeThyraVector::axpby( double alpha, double beta, const VectorOperations &x )
{
  axpby( alpha,
	 beta,
	 *(x.getVectorData()),
	 *getVectorData() );
}


void NativeThyraVector::abs( const VectorOperations &x )
{
    abs( *(x.getVectorData()), *getVectorData() );
}

double NativeThyraVector::min( void ) const
{
  return min( *getVectorData() );
}

double NativeThyraVector::max( void ) const
{
  return max( *getVectorData() );
}

double NativeThyraVector::L1Norm( void ) const
{
  return L1Norm( *getVectorData() );
}


double NativeThyraVector::L2Norm( void ) const
{
  return L2Norm( *getVectorData() );
}


double NativeThyraVector::maxNorm( void ) const
{
  return maxNorm( *getVectorData() );
}


double NativeThyraVector::dot( const VectorOperations &x ) const
{
    return dot( *(x.getVectorData()), *getVectorData() );
}
#endif

void NativeThyraVector::setValuesByLocalID( int num, size_t *indices, const double *vals )
{
    NULL_USE( num );
    NULL_USE( indices );
    NULL_USE( vals );
    AMP_ERROR( "not implemented" );
}


void NativeThyraVector::setLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    NULL_USE( num );
    NULL_USE( indices );
    NULL_USE( vals );
    AMP_ERROR( "not implemented" );
}


void NativeThyraVector::addValuesByLocalID( int num, size_t *indices, const double *vals )
{
    NULL_USE( num );
    NULL_USE( indices );
    NULL_USE( vals );
    AMP_ERROR( "not implemented" );
}


void NativeThyraVector::addLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    NULL_USE( num );
    NULL_USE( indices );
    NULL_USE( vals );
    AMP_ERROR( "not implemented" );
}


void NativeThyraVector::assemble() {}


size_t NativeThyraVector::getLocalSize() const { return d_local; }


size_t NativeThyraVector::getGlobalSize() const { return d_thyraVec->space()->dim(); }


void NativeThyraVector::getLocalValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    NULL_USE( numVals );
    NULL_USE( ndx );
    NULL_USE( vals );
    AMP_ERROR( "not implemented" );
}

void NativeThyraVector::setToScalar( double alpha, VectorData &x )
{
  Thyra::put_scalar<double>( alpha, getThyraVec(x).ptr() );
}

void NativeThyraVector::setRandomValues( VectorData &x )
{
    Thyra::randomize<double>( 0.0, 1.0, getThyraVec(x).ptr() );
}

void NativeThyraVector::setRandomValues( RNG::shared_ptr rng, VectorData &x )
{
    AMP_WARNING("NativeThyraVector::setRandomValue : Not using provided RNG");
    Thyra::randomize<double>( 0.0, 1.0, getThyraVec(x).ptr() );
}

void NativeThyraVector::copy( const VectorData &x, VectorData &y )
{
    Thyra::copy<double>( *( getThyraVec(x) ), getThyraVec(y).ptr() );
}

void NativeThyraVector::scale( double alpha, VectorData &x )
{
  Thyra::scale<double>( alpha, getThyraVec(x).ptr() );
}

void NativeThyraVector::scale( double alpha, const VectorData &x, VectorData &y )
{
    auto src = dynamic_cast<const Vector *>( &x );
    AMP_ASSERT( src != nullptr );
    auto dst = dynamic_cast<Vector *>( &y );
    AMP_ASSERT( dst != nullptr );
    dst->copyVector( src->shared_from_this() );
    Thyra::scale<double>( alpha, getThyraVec(y).ptr() );
}

void NativeThyraVector::add( const VectorData &x, const VectorData &y, VectorData &z )
{
  linearSum( 1.0, x, 1.0, y, z );
}

void NativeThyraVector::subtract( const VectorData &x, const VectorData &y, VectorData &z  )
{
  linearSum( 1.0, x, -1.0, y, z );
}

void NativeThyraVector::multiply( const VectorData &x, const VectorData &y, VectorData &z )
{
    auto xv = getThyraVec( x );
    auto yv = getThyraVec( y );
    auto zv = getThyraVec( z );
    Thyra::put_scalar<double>( 0.0, zv.ptr() );
    Thyra::ele_wise_prod<double>(
        1.0, *xv, *yv, zv.ptr() );
}

void NativeThyraVector::divide( const VectorData &x, const VectorData &y, VectorData &z )
{
    auto xv = getThyraVec( x );
    auto yv = getThyraVec( y );
    auto zv = getThyraVec( z );
    Thyra::put_scalar<double>( 0.0, zv.ptr() );
    Thyra::ele_wise_divide<double>(
        1.0, *xv, *yv, zv.ptr() );
}

void NativeThyraVector::reciprocal( const VectorData &x, VectorData &y )
{
#if TRILINOS_MAJOR_MINOR_VERSION <= 100800
  Thyra::reciprocal<double>( getThyraVec(y).ptr(), *( getThyraVec( x ) ) );
#else
    Thyra::reciprocal<double>( *( getThyraVec( x ) ), getThyraVec(y).ptr() );
#endif
}

void NativeThyraVector::linearSum( double alpha_in,
                                   const VectorData &x,
                                   double beta_in,
                                   const VectorData &y,
				   VectorData &z)
{
    std::vector<double> alpha_vec( 2, 1.0 );
    alpha_vec[0] = alpha_in;
    alpha_vec[1] = beta_in;
    std::vector<Teuchos::Ptr<const Thyra::VectorBase<double>>> vecs( 2 );
    vecs[0] = getThyraVec( x ).ptr();
    vecs[1] = getThyraVec( y ).ptr();
    Teuchos::ArrayView<double> alpha_view( alpha_vec );
    Teuchos::ArrayView<Teuchos::Ptr<const Thyra::VectorBase<double>>> vecs_view( vecs );
    Thyra::linear_combination<double>( alpha_view, vecs_view, 0.0, getThyraVec(z).ptr() );
}

void NativeThyraVector::axpy( double alpha, const VectorData &x, const VectorData &y, VectorData &z )
{
  linearSum( alpha, x, 1.0, y, z );
}

void NativeThyraVector::axpby( double alpha, double beta, const VectorData &x, VectorData &z )
{
  linearSum( alpha, x, beta, z, z );
}

void NativeThyraVector::abs( const VectorData &x, VectorData &y )
{
#if TRILINOS_MAJOR_MINOR_VERSION <= 100800
  Thyra::abs<double>( getThyraVec(y).ptr(), *getThyraVec( x ) );
#else
    Thyra::abs<double>( *getThyraVec( x ), getThyraVec(y).ptr() );
#endif
}

double NativeThyraVector::min( const VectorData &x )  const
{
  return Thyra::min<double>( *getThyraVec( x ) );
}

double NativeThyraVector::max( const VectorData &x )  const
{
  return Thyra::max<double>( *getThyraVec( x ) );
}

double NativeThyraVector::L1Norm( const VectorData &x )  const
{
  return Thyra::norm_1<double>( *getThyraVec( x ) );
}

double NativeThyraVector::L2Norm( const VectorData &x )  const
{
  return Thyra::norm_2<double>( *getThyraVec( x ) );
}

double NativeThyraVector::maxNorm( const VectorData &x )  const
{
  return Thyra::norm_inf<double>( *getThyraVec( x ) );
}

double NativeThyraVector::dot( const VectorData &x, const VectorData &y ) const
{
    return Thyra::dot<double>( *getThyraVec(x), *getThyraVec(y) );
}

#if 0
void NativeThyraVector::addScalar( const VectorData &x, double alpha_in, VectorData &y )
{
}

double NativeThyraVector::localMin( const VectorData &x ) 
{
}

double NativeThyraVector::localMax( const VectorData &x ) 
{
}

double NativeThyraVector::localL1Norm( const VectorData &x ) 
{
}

double NativeThyraVector::localL2Norm( const VectorData &x ) 
{
}

double NativeThyraVector::localMaxNorm( const VectorData &x ) 
{
}

double NativeThyraVector::localDot( const VectorData &x, const VectorData &y )
{
}

double NativeThyraVector::localMinQuotient( const VectorData &x, const VectorData &y )
{
}

double NativeThyraVector::localWrmsNorm( const VectorData &x, const VectorData &y )
{
}

double NativeThyraVector::localWrmsNormMask( const VectorData &x, const VectorData &mask, const VectorData &y )
{
}

bool NativeThyraVector::localEquals( const VectorData &x, const VectorData &y, double tol )
{
}
#endif
} // namespace LinearAlgebra
} // namespace AMP
