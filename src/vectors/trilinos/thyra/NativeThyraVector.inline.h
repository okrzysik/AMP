#include "AMP/vectors/trilinos/thyra/NativeThyraVector.h"

DISABLE_WARNINGS
#include "Thyra_VectorStdOps_def.hpp"
#include "Trilinos_version.h"
ENABLE_WARNINGS


namespace AMP {
namespace LinearAlgebra {

#if 0
inline std::shared_ptr<VectorData> NativeThyraVector::getNewBuffer()
{
    return std::shared_ptr<VectorData>();
}


inline bool NativeThyraVector::sameEngine( VectorEngine &e ) const
{
    return dynamic_cast<NativeThyraVector *>( &e );
}


inline std::shared_ptr<VectorEngine>
    NativeThyraVector::cloneEngine( std::shared_ptr<VectorData> ) const
{
    return std::dynamic_pointer_cast<VectorEngine>( Vector::cloneVector( "engine_clone" ) );
}


inline void NativeThyraVector::swapEngines( std::shared_ptr<VectorEngine> p )
{
    Vector::shared_ptr p2 = std::dynamic_pointer_cast<Vector>( p );
    Vector::swapVectors( p2 );
}

#endif
 
inline void NativeThyraVector::getValuesByLocalID( int numVals, size_t *ndx, double *vals ) const
{
    Vector::getValuesByLocalID( numVals, ndx, vals );
}


inline Vector::shared_ptr NativeThyraVector::getManagedVectorCopy( AMP_MPI comm )
{
    NULL_USE( comm );
    AMP_ERROR( "Not programmed yet" );
    return Vector::shared_ptr();
}

#if 0
inline AMP_MPI NativeThyraVector::getComm() const { return Vector::getComm(); }
#endif
 
inline Vector::shared_ptr NativeThyraVector::getManagedVectorDuplicate( AMP_MPI comm )
{
    NULL_USE( comm );

    AMP_ERROR( "Not programmed yet" );
    return Vector::shared_ptr();
}


inline void NativeThyraVector::aliasVector( Vector & ) { AMP_ERROR( "not implemented" ); }


inline void NativeThyraVector::swapVectors( Vector & ) { AMP_ERROR( "not implemented" ); }


inline Teuchos::RCP<const Thyra::VectorBase<double>>
NativeThyraVector::getThyraVec( const VectorOperations &v )
{
    auto vec = dynamic_cast<const Vector *>( &v );
    AMP_ASSERT( vec != nullptr );
    auto vec2 = std::dynamic_pointer_cast<const ThyraVector>(
        ThyraVector::constView( vec->shared_from_this() ) );
    AMP_ASSERT( vec2 != nullptr );
    return vec2->getVec();
}


inline Teuchos::RCP<const Thyra::VectorBase<double>>
NativeThyraVector::getThyraVec( const Vector::const_shared_ptr &vec )
{
    auto vec2 = std::dynamic_pointer_cast<const ThyraVector>( ThyraVector::constView( vec ) );
    AMP_ASSERT( vec2 != nullptr );
    return vec2->getVec();
}


inline void NativeThyraVector::copy( const VectorOperations &src_vec )
{
    Thyra::copy<double>( *( getThyraVec( src_vec ) ), d_thyraVec.ptr() );
}


inline void NativeThyraVector::setToScalar( double alpha )
{
    Thyra::put_scalar<double>( alpha, d_thyraVec.ptr() );
}


inline void NativeThyraVector::scale( double alpha, const VectorOperations &x )
{
    auto vec = dynamic_cast<const Vector *>( &x );
    AMP_ASSERT( vec != nullptr );
    copyVector( vec->shared_from_this() );
    Thyra::scale<double>( alpha, d_thyraVec.ptr() );
}


inline void NativeThyraVector::scale( double alpha )
{
    Thyra::scale<double>( alpha, d_thyraVec.ptr() );
}


inline void NativeThyraVector::add( const VectorOperations &x, const VectorOperations &y )
{
    linearSum( 1.0, x, 1.0, y );
}


inline void NativeThyraVector::subtract( const VectorOperations &x, const VectorOperations &y )
{
    linearSum( 1.0, x, -1.0, y );
}


inline void NativeThyraVector::multiply( const VectorOperations &x, const VectorOperations &y )
{
    Thyra::put_scalar<double>( 0.0, d_thyraVec.ptr() );
    Thyra::ele_wise_prod<double>(
        1.0, *( getThyraVec( x ) ), *( getThyraVec( y ) ), d_thyraVec.ptr() );
}


inline void NativeThyraVector::divide( const VectorOperations &x, const VectorOperations &y )
{
    Thyra::put_scalar<double>( 0.0, d_thyraVec.ptr() );
    Thyra::ele_wise_divide<double>(
        1.0, *( getThyraVec( x ) ), *( getThyraVec( y ) ), d_thyraVec.ptr() );
}


inline void NativeThyraVector::reciprocal( const VectorOperations &x )
{
#if TRILINOS_MAJOR_MINOR_VERSION <= 100800
    Thyra::reciprocal<double>( d_thyraVec.ptr(), *( getThyraVec( x ) ) );
#else
    Thyra::reciprocal<double>( *( getThyraVec( x ) ), d_thyraVec.ptr() );
#endif
}


inline void NativeThyraVector::linearSum( double alpha,
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


inline void
NativeThyraVector::axpy( double alpha, const VectorOperations &x, const VectorOperations &y )
{
    linearSum( alpha, x, 1.0, y );
}


inline void NativeThyraVector::axpby( double alpha, double beta, const VectorOperations &x )
{
    linearSum( alpha, x, beta, *this );
}


inline void NativeThyraVector::abs( const VectorOperations &x )
{
#if TRILINOS_MAJOR_MINOR_VERSION <= 100800
    Thyra::abs<double>( d_thyraVec.ptr(), *getThyraVec( x ) );
#else
    Thyra::abs<double>( *getThyraVec( x ), d_thyraVec.ptr() );
#endif
}


inline double NativeThyraVector::min( void ) const { return Thyra::min<double>( *d_thyraVec ); }


inline double NativeThyraVector::max( void ) const { return Thyra::max<double>( *d_thyraVec ); }


inline void NativeThyraVector::setRandomValues( void )
{
    Thyra::randomize<double>( 0.0, 1.0, d_thyraVec.ptr() );
}


inline double NativeThyraVector::L1Norm( void ) const
{
    return Thyra::norm_1<double>( *d_thyraVec );
}


inline double NativeThyraVector::L2Norm( void ) const
{
    return Thyra::norm_2<double>( *d_thyraVec );
}


inline double NativeThyraVector::maxNorm( void ) const
{
    return Thyra::norm_inf<double>( *d_thyraVec );
}


inline double NativeThyraVector::dot( const VectorOperations &x ) const
{
    return Thyra::dot<double>( *getThyraVec( x ), *d_thyraVec );
}


inline void NativeThyraVector::setValuesByLocalID( int num, size_t *indices, const double *vals )
{
    NULL_USE( num );
    NULL_USE( indices );
    NULL_USE( vals );
    AMP_ERROR( "not implemented" );
}


inline void
NativeThyraVector::setLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    NULL_USE( num );
    NULL_USE( indices );
    NULL_USE( vals );
    AMP_ERROR( "not implemented" );
}


inline void NativeThyraVector::addValuesByLocalID( int num, size_t *indices, const double *vals )
{
    NULL_USE( num );
    NULL_USE( indices );
    NULL_USE( vals );
    AMP_ERROR( "not implemented" );
}


inline void
NativeThyraVector::addLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    NULL_USE( num );
    NULL_USE( indices );
    NULL_USE( vals );
    AMP_ERROR( "not implemented" );
}


inline void NativeThyraVector::assemble() {}


inline size_t NativeThyraVector::getLocalSize() const { return d_local; }


inline size_t NativeThyraVector::getGlobalSize() const { return d_thyraVec->space()->dim(); }


inline void
NativeThyraVector::getLocalValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    NULL_USE( numVals );
    NULL_USE( ndx );
    NULL_USE( vals );
    AMP_ERROR( "not implemented" );
}


} // namespace LinearAlgebra
} // namespace AMP
