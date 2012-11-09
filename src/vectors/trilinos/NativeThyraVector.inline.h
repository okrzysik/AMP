#include "vectors/trilinos/NativeThyraVector.h"


#include "Trilinos_version.h"
#include "Thyra_VectorStdOps_def.hpp"


namespace AMP {
namespace LinearAlgebra {


inline VectorEngine::BufferPtr   NativeThyraVector::getNewBuffer ()
{
    return VectorEngine::BufferPtr ();
}


inline bool NativeThyraVector::sameEngine ( VectorEngine &e ) const
{
    return e.isA<NativeThyraVector> ();
}


inline VectorEngine::shared_ptr  NativeThyraVector::cloneEngine ( BufferPtr ) const
{
    return boost::dynamic_pointer_cast<VectorEngine> ( Vector::cloneVector ( "engine_clone" ) );
}


inline void NativeThyraVector::swapEngines ( VectorEngine::shared_ptr p )
{
    Vector::shared_ptr p2 =boost::dynamic_pointer_cast<Vector> ( p );
    Vector::swapVectors ( p2 );
}


inline void *NativeThyraVector::getDataBlock ( size_t i )
{
    return static_cast<void *>(getRawDataBlock<double> ( i ));
}


inline const void *NativeThyraVector::getDataBlock ( size_t i ) const
{
    return static_cast<const void *>(getRawDataBlock<double> ( i ));
}


inline void NativeThyraVector::getValuesByLocalID ( int numVals , size_t *ndx , double *vals ) const
{
    INCREMENT_COUNT("Virtual");
    Vector::getValuesByLocalID ( numVals , ndx , vals );
}


inline void NativeThyraVector::copyOutRawData ( double **out )
{
    *out = new double [getLocalSize()];
    std::copy ( getRawDataBlock<double> ( 0 ) ,
                getRawDataBlock<double> ( 0 ) + getLocalSize() ,
                *out );
}


inline Vector::shared_ptr  NativeThyraVector::getManagedVectorCopy ( AMP_MPI comm )
{
    AMP_ERROR("Not programmed yet");
    return Vector::shared_ptr();
}


inline AMP_MPI NativeThyraVector::getComm () const
{
    return Vector::getComm();
}


inline Vector::shared_ptr  NativeThyraVector::getManagedVectorDuplicate ( AMP_MPI  comm )
{
    AMP_ERROR("Not programmed yet");
    return Vector::shared_ptr();
}


inline void NativeThyraVector::aliasVector(Vector &)
{
    AMP_ERROR( "not implemented" );
}


inline void NativeThyraVector::swapVectors(Vector &other)
{
    AMP_ERROR( "not implemented" );
}


inline Teuchos::RCP<const Thyra::VectorBase<double> >  NativeThyraVector::getThyraVec( const VectorOperations &v )
{
    boost::shared_ptr<const ThyraVector> v2 = boost::dynamic_pointer_cast<const ThyraVector>( 
        ThyraVector::constView( v.castTo<const Vector>().shared_from_this() ) );
    AMP_ASSERT(v2!=NULL);
    return v2->getVec();
}


inline Teuchos::RCP<const Thyra::VectorBase<double> >  NativeThyraVector::getThyraVec( const Vector::const_shared_ptr &v )
{
    boost::shared_ptr<const ThyraVector> v2 = 
        boost::dynamic_pointer_cast<const ThyraVector>( ThyraVector::constView( v ) );
    AMP_ASSERT(v2!=NULL);
    return v2->getVec();
}


inline void NativeThyraVector::copyVector(const Vector::const_shared_ptr &src_vec)
{
    Thyra::copy<double>( *(getThyraVec(src_vec)), d_thyraVec );
}


inline void  NativeThyraVector::setToScalar(double alpha)
{
    Thyra::put_scalar<double>( alpha, d_thyraVec );
}


inline void  NativeThyraVector::scale(double alpha, const VectorOperations &x)
{
    copyVector( x.castTo<const Vector>().shared_from_this() );
    Thyra::scale<double>( alpha, d_thyraVec );
}


inline void  NativeThyraVector::scale(double alpha)
{
    Thyra::scale<double>( alpha, d_thyraVec );
}


inline void  NativeThyraVector::add(const VectorOperations &x, const VectorOperations &y)
{
    linearSum( 1.0, x, 1.0, y );
}


inline void  NativeThyraVector::subtract(const VectorOperations &x, const VectorOperations &y)
{
    linearSum( 1.0, x, -1.0, y );
}


inline void  NativeThyraVector::multiply( const VectorOperations &x, const VectorOperations &y)
{
    Thyra::put_scalar<double>( 0.0, d_thyraVec );
    Thyra::ele_wise_prod<double>( 1.0, *(getThyraVec(x)), *(getThyraVec(y)), d_thyraVec );
}


inline void  NativeThyraVector::divide( const VectorOperations &x, const VectorOperations &y)
{
    Thyra::put_scalar<double>( 0.0, d_thyraVec );
    Thyra::ele_wise_divide<double>( 1.0, *(getThyraVec(x)), *(getThyraVec(y)), d_thyraVec );
}


inline void  NativeThyraVector::reciprocal(const VectorOperations &x)
{
    #if TRILINOS_MAJOR_MINOR_VERSION <= 100800
        Thyra::reciprocal<double>( d_thyraVec, *(getThyraVec(x)) );
    #else
        Thyra::reciprocal<double>( *(getThyraVec(x)), d_thyraVec );
    #endif
}


inline void  NativeThyraVector::linearSum(double alpha, const VectorOperations &x,
       double beta, const VectorOperations &y)
{
    std::vector<double> alpha_vec(2,1.0);
    alpha_vec[0] = alpha;
    alpha_vec[1] = beta;
    std::vector<Teuchos::Ptr<const Thyra::VectorBase<double> > > vecs(2);
    vecs[0] = getThyraVec(x);
    vecs[1] = getThyraVec(y);
    Teuchos::ArrayView<double> alpha_view( alpha_vec );
    Teuchos::ArrayView<Teuchos::Ptr<const Thyra::VectorBase<double> > > vecs_view( vecs );
    Thyra::linear_combination<double>(	alpha_view, vecs_view, 0.0, d_thyraVec );
}


inline void  NativeThyraVector::axpy(double alpha, const VectorOperations &x, const VectorOperations &y)
{
    linearSum( alpha, x, 1.0, y );
}


inline void  NativeThyraVector::axpby(double alpha, double beta, const VectorOperations &x)
{
    linearSum( alpha, x, beta, *this );
}


inline void  NativeThyraVector::abs(const VectorOperations &x)
{
    #if TRILINOS_MAJOR_MINOR_VERSION <= 100800
       Thyra::abs<double>( d_thyraVec, *getThyraVec(x) );
    #else
        Thyra::abs<double>( *getThyraVec(x), d_thyraVec );
    #endif
}


inline double NativeThyraVector::min(void) const
{
    return Thyra::min<double>( *d_thyraVec );
}


inline double NativeThyraVector::max(void) const
{
    return Thyra::max<double>( *d_thyraVec );
}


inline void  NativeThyraVector::setRandomValues(void)
{
    Thyra::randomize<double>( 0.0, 1.0, d_thyraVec );
}


inline double NativeThyraVector::L1Norm(void) const
{
    return Thyra::norm_1<double>( *d_thyraVec );
}


inline double NativeThyraVector::L2Norm(void) const
{
    return Thyra::norm_2<double>( *d_thyraVec );
}


inline double NativeThyraVector::maxNorm(void) const
{
    return Thyra::norm_inf<double>( *d_thyraVec );
}


inline double NativeThyraVector::dot(const VectorOperations &x) const
{
    return Thyra::dot<double>( *getThyraVec(x), *d_thyraVec );
}


inline void  NativeThyraVector::setValuesByLocalID(int num , size_t *indices , const double *vals)
{
    AMP_ERROR( "not implemented" );
}


inline void  NativeThyraVector::setLocalValuesByGlobalID(int num, size_t *indices , const double *vals)
{
    AMP_ERROR( "not implemented" );
}


inline void  NativeThyraVector::addValuesByLocalID(int num , size_t *indices , const double *vals)
{
    AMP_ERROR( "not implemented" );
}


inline void  NativeThyraVector::addLocalValuesByGlobalID(int num, size_t *indices , const double *vals)
{
    AMP_ERROR( "not implemented" );
}


inline void  NativeThyraVector::assemble()
{
}


inline size_t NativeThyraVector::getLocalSize() const
{
    return d_local;
}


inline size_t NativeThyraVector::getGlobalSize() const
{
    return d_thyraVec->space()->dim();
}


inline void NativeThyraVector::getLocalValuesByGlobalID ( int numVals , size_t *ndx , double *vals ) const
{
    AMP_ERROR( "not implemented" );
}


}
}

