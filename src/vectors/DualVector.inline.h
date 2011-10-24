
#include "math.h"

/// \cond UNDOCUMENTED

namespace AMP {
namespace LinearAlgebra {

  inline
  void  DualVector::selectInto ( const VectorSelector &s , Vector::shared_ptr p ) 
  {
    d_pVector1->selectInto ( s , p );
    d_pVector2->selectInto ( s , p );
  }

  inline
  AMP_MPI  DualVector::getComm () const 
  { 
    return d_pVector1->getComm(); 
  }

  inline
  DualVector::~DualVector ()
  {
  }

  inline
  Vector::shared_ptr  DualVector::create ( Vector::shared_ptr v1 , Vector::shared_ptr v2 , const std::string &names )
  {
    Variable::shared_ptr n ( new DualVariable ( v1->getVariable() , v2->getVariable() , names ) );
    return Vector::shared_ptr ( new DualVector ( v1 , v2 , n ) );
  }

  inline
  Vector::shared_ptr  DualVector::create ( Vector::shared_ptr v1 , Vector::shared_ptr v2 , Variable::shared_ptr names )
  {
    return Vector::shared_ptr ( new DualVector ( v1 , v2 , names ) );
  }


  inline
  size_t DualVector::numberOfDataBlocks () const 
  { 
    size_t a = d_pVector1->numberOfDataBlocks ();
    size_t b = d_pVector2->numberOfDataBlocks ();
    return a + b; 
  }

  inline
  size_t DualVector::sizeOfDataBlock ( size_t i ) const
  {
    size_t retVal;
    if ( i >= d_pVector1->numberOfDataBlocks() )
      retVal = d_pVector2->sizeOfDataBlock ( i - d_pVector1->numberOfDataBlocks() );
    else
      retVal = d_pVector1->sizeOfDataBlock ( i );
    return retVal;
  }

  inline
  void * DualVector::getDataBlock ( size_t i ) 
  { 
    return getRawDataBlockAsVoid ( i ); 
  }

  inline
  const void * DualVector::getDataBlock ( size_t i ) const 
  { 
    return getRawDataBlockAsVoid ( i ); 
  }

  inline
  VectorEngine::BufferPtr  DualVector::getNewBuffer () 
  { 
    return BufferPtr (); 
  }

  inline
  bool  DualVector::sameEngine ( VectorEngine &p ) const 
  { 
    return p.isA<DualVector>(); 
  }

  inline
  void  DualVector::swapEngines ( VectorEngine::shared_ptr p ) 
  { 
    swapVectors ( p->castTo<DualVector>() ); 
  }

  inline
  VectorEngine::shared_ptr  DualVector::cloneEngine ( BufferPtr ) const 
  { 
    return boost::dynamic_pointer_cast<VectorEngine> ( Vector::cloneVector( "engine_var" ) ); 
  }

  inline
  const Vector::shared_ptr  &DualVector::getVector1 ( const VectorOperations &rhs ) 
  {
    return rhs.castTo<DualVector>().d_pVector1;
  }
  
  inline
  const Vector::shared_ptr  &DualVector::getVector2 ( const VectorOperations &rhs )  
  {
    return rhs.castTo<DualVector>().d_pVector2;
  }

  inline
  Vector::shared_ptr  &DualVector::getVector1 ( VectorOperations &rhs )  
  {
    return rhs.castTo<DualVector>().d_pVector1;
  }
  
  inline
  Vector::shared_ptr  &DualVector::getVector2 ( VectorOperations &rhs ) 
  {
    return rhs.castTo<DualVector>().d_pVector2;
  }

  inline
  void * DualVector::getRawDataBlockAsVoid ( size_t i )
  {
    void * retVal;
    if ( i >= d_pVector1->numberOfDataBlocks() )
      retVal = static_cast<void *> ( d_pVector2->getRawDataBlock<double>( i - d_pVector1->numberOfDataBlocks() ) );
    else
      retVal = static_cast<void *> ( d_pVector1->getRawDataBlock<double> ( i ) );
    return retVal;
  }

  inline
  const void * DualVector::getRawDataBlockAsVoid ( size_t i ) const
  {
    void * retVal;
    if ( i >= d_pVector1->numberOfDataBlocks() )
      retVal = static_cast<void *> ( d_pVector2->getRawDataBlock<double>( i - d_pVector1->numberOfDataBlocks() ) );
    else
      retVal = static_cast<void *> ( d_pVector1->getRawDataBlock<double> ( i ) );
    return retVal;
  }

  inline
  void DualVector::copyVector ( const Vector &src )
  {
     d_pVector1->copyVector ( getVector1( src ) );
     d_pVector2->copyVector ( getVector2( src ) );
  }

  inline
  void DualVector::swapVectors(Vector &other)
  {
    d_pVector1->swapVectors ( getVector1( other ) );
    d_pVector2->swapVectors ( getVector2( other ) );
  }


  inline
  Vector::shared_ptr DualVector::cloneVector(const Variable::shared_ptr name) const
  {
    return Vector::shared_ptr ( new DualVector ( d_pVector1->cloneVector ( d_pVector1->getVariable() ) , 
                                         d_pVector2->cloneVector ( d_pVector2->getVariable() ) , name ) );
  }

  inline
  void DualVector::aliasVector (Vector &other)
  {
    d_pVector1->aliasVector ( getVector1 ( other ) );
    d_pVector2->aliasVector ( getVector2 ( other ) );
  }

  inline
  void DualVector::setToScalar ( double alpha )
  {
    d_pVector1->setToScalar ( alpha );
    d_pVector2->setToScalar ( alpha );
  }

  inline
  void DualVector::scale ( double alpha , const VectorOperations &x )
  {
    d_pVector1->scale ( alpha , getVector1 ( x ) );
    d_pVector2->scale ( alpha , getVector2 ( x ) );
  }

  inline
  void DualVector::scale ( double alpha )
  {
    d_pVector1->scale ( alpha );
    d_pVector2->scale ( alpha );
  }

  inline
  void DualVector::add ( const VectorOperations &x , const VectorOperations &y )
  {
    d_pVector1->add ( getVector1 ( x ) , getVector1 ( y ) );
    d_pVector2->add ( getVector2 ( x ) , getVector2 ( y ) );
  }

  inline
  void DualVector::subtract ( const VectorOperations &x , const VectorOperations &y )
  {
    d_pVector1->subtract ( getVector1 ( x ) , getVector1 ( y ) );
    d_pVector2->subtract ( getVector2 ( x ) , getVector2 ( y ) );
  }

  inline
  void DualVector::multiply ( const VectorOperations &x , const VectorOperations &y )
  {
    d_pVector1->multiply ( getVector1 ( x ) , getVector1 ( y ) );
    d_pVector2->multiply ( getVector2 ( x ) , getVector2 ( y ) );
  }

  inline
  void DualVector::divide ( const VectorOperations &x , const VectorOperations &y )
  {
    d_pVector1->divide ( getVector1 ( x ) , getVector1 ( y ) );
    d_pVector2->divide ( getVector2 ( x ) , getVector2 ( y ) );
  }

  inline
  void DualVector::reciprocal ( const VectorOperations &x )
  {
    d_pVector1->reciprocal ( getVector1 ( x ) );
    d_pVector2->reciprocal ( getVector2 ( x ) );
  }

  inline
  void DualVector::linearSum(double alpha, const VectorOperations &x, double beta, const VectorOperations &y)
  {
    d_pVector1->linearSum ( alpha , getVector1 ( x ) , beta , getVector1 ( y ) );
    d_pVector2->linearSum ( alpha , getVector2 ( x ) , beta , getVector2 ( y ) );
  }

  inline
  void DualVector::axpy(double alpha, const VectorOperations &x, const VectorOperations &y)
  {
    d_pVector1->axpy ( alpha , getVector1 ( x ) , getVector1 ( y ) );
    d_pVector2->axpy ( alpha , getVector2 ( x ) , getVector2 ( y ) );
  }

  inline
  void DualVector::axpby(double alpha, double beta, const VectorOperations &x)
  {
    d_pVector1->axpby ( alpha , beta , getVector1 ( x ) );
    d_pVector2->axpby ( alpha , beta , getVector2 ( x ) );
  }

  inline
  void DualVector::abs ( const VectorOperations &x )
  {
    d_pVector1->abs ( getVector1 ( x ) );
    d_pVector2->abs ( getVector2 ( x ) );
  }

  inline
  double DualVector::min(void) const
  {
    return std::min ( d_pVector1->min() , d_pVector2->min() );
  }

  inline
  double DualVector::max(void) const
  {
    return std::max ( d_pVector1->max() , d_pVector2->max() );
  }

  inline
  void DualVector::setRandomValues ()
  {
    d_pVector1->setRandomValues ();
    d_pVector2->setRandomValues ();
  }

  inline
  void DualVector::makeConsistent ( ScatterType t )
  {
    d_pVector1->makeConsistent ( t );
    d_pVector2->makeConsistent ( t );
  }

  inline
  void DualVector::assemble () 
  {
    d_pVector1->assemble ( );
    d_pVector2->assemble ( );
  }

  inline
  double DualVector::L1Norm () const
  {
    return d_pVector1->L1Norm() + d_pVector2->L1Norm();
  }

  inline
  double DualVector::L2Norm () const
  {
    double n1 = d_pVector1->L2Norm();
    double n2 = d_pVector2->L2Norm();
    return sqrt ( n1 * n1 + n2 * n2 );
  }

  inline
  double DualVector::maxNorm () const
  {
    return std::max ( d_pVector1->maxNorm() , d_pVector2->maxNorm() );
  }

  inline
  double DualVector::dot ( const VectorOperations &rhs ) const
  {
    return d_pVector1->dot ( getVector1 ( rhs ) ) + d_pVector2->dot ( getVector2 ( rhs ) );
  }

  inline
  size_t DualVector::getLocalSize () const
  {
    return d_pVector1->getLocalSize() + d_pVector2->getLocalSize();
  }

  inline
  size_t DualVector::getGlobalSize () const
  {
    return d_pVector1->getGlobalSize() + d_pVector2->getGlobalSize();
  }

  inline
  size_t DualVector::getGhostSize () const
  {
    return d_pVector1->getGhostSize() + d_pVector2->getGhostSize();
  }

  inline
  void  DualVector::putRawData ( double *in )
  {
    d_pVector1->putRawData ( in );
    d_pVector2->putRawData ( in + d_pVector1->getLocalSize() );
  }

  inline
  void  DualVector::copyOutRawData ( double **out )
  {
    *out = new double [ getLocalSize() ];
    d_pVector1->copyOutRawData ( out );
    d_pVector2->copyOutRawData ( out + d_pVector2->getLocalSize() );
  }

}
}

/// \endcond
