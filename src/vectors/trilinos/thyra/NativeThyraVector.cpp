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
    : NativeVector(), ThyraVector()
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


Vector::shared_ptr NativeThyraVector::getManagedVectorCopy( AMP_MPI comm )
{
    NULL_USE( comm );
    AMP_ERROR( "Not programmed yet" );
    return Vector::shared_ptr();
}

Vector::shared_ptr NativeThyraVector::getManagedVectorDuplicate( AMP_MPI comm )
{
    NULL_USE( comm );

    AMP_ERROR( "Not programmed yet" );
    return Vector::shared_ptr();
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


void NativeThyraVector::copy( const VectorOperations &src_vec )
{
    Thyra::copy<double>( *( getThyraVec( src_vec ) ), d_thyraVec.ptr() );
}


void NativeThyraVector::setToScalar( double alpha )
{
    Thyra::put_scalar<double>( alpha, d_thyraVec.ptr() );
}


void NativeThyraVector::scale( double alpha, const VectorOperations &x )
{
    auto vec = dynamic_cast<const Vector *>( &x );
    AMP_ASSERT( vec != nullptr );
    copyVector( vec->shared_from_this() );
    Thyra::scale<double>( alpha, d_thyraVec.ptr() );
}


void NativeThyraVector::scale( double alpha ) { Thyra::scale<double>( alpha, d_thyraVec.ptr() ); }


void NativeThyraVector::add( const VectorOperations &x, const VectorOperations &y )
{
    linearSum( 1.0, x, 1.0, y );
}


void NativeThyraVector::subtract( const VectorOperations &x, const VectorOperations &y )
{
    linearSum( 1.0, x, -1.0, y );
}


void NativeThyraVector::multiply( const VectorOperations &x, const VectorOperations &y )
{
    Thyra::put_scalar<double>( 0.0, d_thyraVec.ptr() );
    Thyra::ele_wise_prod<double>(
        1.0, *( getThyraVec( x ) ), *( getThyraVec( y ) ), d_thyraVec.ptr() );
}


void NativeThyraVector::divide( const VectorOperations &x, const VectorOperations &y )
{
    Thyra::put_scalar<double>( 0.0, d_thyraVec.ptr() );
    Thyra::ele_wise_divide<double>(
        1.0, *( getThyraVec( x ) ), *( getThyraVec( y ) ), d_thyraVec.ptr() );
}


void NativeThyraVector::reciprocal( const VectorOperations &x )
{
#if TRILINOS_MAJOR_MINOR_VERSION <= 100800
    Thyra::reciprocal<double>( d_thyraVec.ptr(), *( getThyraVec( x ) ) );
#else
    Thyra::reciprocal<double>( *( getThyraVec( x ) ), d_thyraVec.ptr() );
#endif
}


void NativeThyraVector::linearSum( double alpha,
                                   const VectorOperations &x,
                                   double beta,
                                   const VectorOperations &y )
{
    std::vector<double> alpha_vec( 2, 1.0 );
    alpha_vec[0] = alpha;
    alpha_vec[1] = beta;
    std::vector<Teuchos::Ptr<const Thyra::VectorBase<double>>> vecs( 2 );
    vecs[0] = getThyraVec( x ).ptr();
    vecs[1] = getThyraVec( y ).ptr();
    Teuchos::ArrayView<double> alpha_view( alpha_vec );
    Teuchos::ArrayView<Teuchos::Ptr<const Thyra::VectorBase<double>>> vecs_view( vecs );
    Thyra::linear_combination<double>( alpha_view, vecs_view, 0.0, d_thyraVec.ptr() );
}


void NativeThyraVector::axpy( double alpha, const VectorOperations &x, const VectorOperations &y )
{
    linearSum( alpha, x, 1.0, y );
}


void NativeThyraVector::axpby( double alpha, double beta, const VectorOperations &x )
{
    linearSum( alpha, x, beta, *this );
}


void NativeThyraVector::abs( const VectorOperations &x )
{
#if TRILINOS_MAJOR_MINOR_VERSION <= 100800
    Thyra::abs<double>( d_thyraVec.ptr(), *getThyraVec( x ) );
#else
    Thyra::abs<double>( *getThyraVec( x ), d_thyraVec.ptr() );
#endif
}


double NativeThyraVector::min( void ) const { return Thyra::min<double>( *d_thyraVec ); }


double NativeThyraVector::max( void ) const { return Thyra::max<double>( *d_thyraVec ); }


void NativeThyraVector::setRandomValues( void )
{
    Thyra::randomize<double>( 0.0, 1.0, d_thyraVec.ptr() );
}


double NativeThyraVector::L1Norm( void ) const { return Thyra::norm_1<double>( *d_thyraVec ); }


double NativeThyraVector::L2Norm( void ) const { return Thyra::norm_2<double>( *d_thyraVec ); }


double NativeThyraVector::maxNorm( void ) const { return Thyra::norm_inf<double>( *d_thyraVec ); }


double NativeThyraVector::dot( const VectorOperations &x ) const
{
    return Thyra::dot<double>( *getThyraVec( x ), *d_thyraVec );
}


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

} // namespace LinearAlgebra
} // namespace AMP
