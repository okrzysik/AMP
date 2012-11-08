#include "vectors/trilinos/NativeThyraVector.h"


namespace AMP {
namespace LinearAlgebra {
  

inline NativeThyraVectorParameters::NativeThyraVectorParameters( Teuchos::RCP<Thyra::VectorSpaceBase<double> > v )
{
    d_InVec = v;
}


inline size_t NativeThyraVector::numberOfDataBlocks () const 
{ 
    return 1; 
}


inline size_t NativeThyraVector::sizeOfDataBlock ( size_t i ) const
{
    if ( i != 0 )
      return 0;
    return getLocalSize();
}


inline VectorEngine::BufferPtr   NativeThyraVector::getNewBuffer ()
{
    return VectorEngine::BufferPtr ();
}


inline bool        NativeThyraVector::sameEngine ( VectorEngine &e ) const
{
    return e.isA<NativeThyraVector> ();
}


inline VectorEngine::shared_ptr  NativeThyraVector::cloneEngine ( BufferPtr ) const
{
    return boost::dynamic_pointer_cast<VectorEngine> ( Vector::cloneVector ( "engine_clone" ) );
}


inline void        NativeThyraVector::swapEngines ( VectorEngine::shared_ptr p )
{
    Vector::shared_ptr p2 =boost::dynamic_pointer_cast<Vector> ( p );
    Vector::swapVectors ( p2 );
}


inline void       *NativeThyraVector::getDataBlock ( size_t i )
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


inline void NativeThyraVector::copyVector(const Vector::const_shared_ptr &src_vec)
{
    AMP_ERROR( "not implemented" );
}


inline void  NativeThyraVector::setToScalar(double alpha)
{
    AMP_ERROR( "not implemented" );
}


inline void  NativeThyraVector::scale(double alpha, const VectorOperations &x)
{
    AMP_ERROR( "not implemented" );
}


inline void  NativeThyraVector::scale(double alpha)
{
    AMP_ERROR( "not implemented" );
}


inline void  NativeThyraVector::add(const VectorOperations &x, const VectorOperations &y)
{
    AMP_ERROR( "not implemented" );
}


inline void  NativeThyraVector::subtract(const VectorOperations &x, const VectorOperations &y)
{
    AMP_ERROR( "not implemented" );

}


inline void  NativeThyraVector::multiply( const VectorOperations &x, const VectorOperations &y)
{
    AMP_ERROR( "not implemented" );
}


inline void  NativeThyraVector::divide( const VectorOperations &x, const VectorOperations &y)
{
    AMP_ERROR( "not implemented" );
}


inline void  NativeThyraVector::reciprocal(const VectorOperations &x)
{
    AMP_ERROR( "not implemented" );
}


inline void  NativeThyraVector::linearSum(double alpha, const VectorOperations &x,
       double beta, const VectorOperations &y)
{
    AMP_ERROR( "not implemented" );
}


inline void  NativeThyraVector::axpy(double alpha, const VectorOperations &x, const VectorOperations &y)
{
    AMP_ERROR( "not implemented" );
}


inline void  NativeThyraVector::axpby(double alpha, double beta, const VectorOperations &x)
{
    AMP_ERROR( "not implemented" );
}


inline void  NativeThyraVector::abs(const VectorOperations &x)
{
    AMP_ERROR( "not implemented" );
}


inline double NativeThyraVector::min(void) const
{
    AMP_ERROR( "not implemented" );
    return 0;
}


inline double NativeThyraVector::max(void) const
{
    AMP_ERROR( "not implemented" );
    return 0;
}


inline void  NativeThyraVector::setRandomValues(void)
{
    AMP_ERROR( "not implemented" );
}


inline double NativeThyraVector::L1Norm(void) const
{
    AMP_ERROR( "not implemented" );
    return 0;
}


inline double NativeThyraVector::L2Norm(void) const
{
    AMP_ERROR( "not implemented" );
    return 0;
}


inline double NativeThyraVector::maxNorm(void) const
{
    AMP_ERROR( "not implemented" );
    return 0;
}


inline double NativeThyraVector::dot(const VectorOperations &x) const
{
    AMP_ERROR( "not implemented" );
    return 0;
}


inline double NativeThyraVector::localL1Norm(void) const
{
    AMP_ERROR( "not implemented" );
    return 0;
}


inline double NativeThyraVector::localL2Norm(void) const
{
    AMP_ERROR( "not implemented" );
    return 0;
}


inline double NativeThyraVector::localMaxNorm(void) const
{
    AMP_ERROR( "not implemented" );
    return 0;
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
    AMP_ERROR( "not implemented" );
    return 0;
}


inline size_t NativeThyraVector::getGlobalSize() const
{
    AMP_ERROR( "not implemented" );
    return 0;
}


inline void NativeThyraVector::getLocalValuesByGlobalID ( int numVals , size_t *ndx , double *vals ) const
{
    AMP_ERROR( "not implemented" );
}


}
}

