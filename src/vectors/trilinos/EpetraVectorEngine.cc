#include "utils/Utilities.h"

#include "EpetraVectorEngine.h"

#ifdef USE_MPI
    #include <Epetra_MpiComm.h>
#else
    #include <Epetra_SerialComm.h>
#endif

namespace AMP {
namespace LinearAlgebra {

  Epetra_Map  &EpetraVectorEngineParameters::getEpetraMap()
  {
    #ifdef USE_MPI
        Epetra_MpiComm  comm = d_comm.getCommunicator();
    #else
        Epetra_SerialComm  comm;
    #endif
    if ( d_emap.get() == 0 )
      d_emap = boost::shared_ptr<Epetra_Map> ( new Epetra_Map ( getGlobalSize() , getLocalSize() , &*begin() , 0 , comm ) );
    return *d_emap;
  }

  EpetraVectorEngine::EpetraVectorEngine ( VectorEngineParameters::shared_ptr  alias , BufferPtr buf )
    : d_epetraVector ( View , alias->castTo<EpetraVectorEngineParameters>().getEpetraMap() , &*(buf->begin()) )
  {
    d_Params = alias;
  }

  VectorEngine::BufferPtr  EpetraVectorEngine::getNewBuffer ()
  { 
    BufferPtr retval ( new std::vector<double> ( getLocalSize() ) ); 
//    d_epetraVector.ResetView ( &*(retval->begin()) );
    return retval;
  }

  void EpetraVectorEngine::swapEngines ( VectorEngine::shared_ptr x )
  {
    double * my_pointer; 
    double * oth_pointer; 
    getEpetra_Vector().ExtractView ( &my_pointer );
    x->castTo<EpetraVectorEngine>().getEpetra_Vector().ExtractView ( &oth_pointer );
    x->castTo<EpetraVectorEngine>().getEpetra_Vector().ResetView ( my_pointer );
    getEpetra_Vector().ResetView ( oth_pointer );
  }

  void EpetraVectorEngine::setToScalar(const double alpha)
  {
    getEpetra_Vector().PutScalar ( alpha );
  }

  void EpetraVectorEngine::scale(double alpha, const VectorOperations &x) 
  {
    getEpetra_Vector().Scale ( alpha , x.castTo<EpetraVectorEngine>().getEpetra_Vector() );
  }

  void EpetraVectorEngine::scale(double alpha)
  {
    getEpetra_Vector().Scale ( alpha );
  }

  VectorEngine::shared_ptr EpetraVectorEngine::cloneEngine ( BufferPtr p ) const
  {
    return shared_ptr ( new EpetraVectorEngine ( d_Params , p ) );
  }

  void EpetraVectorEngine::add(const VectorOperations &x, const VectorOperations &y)
  {
    getEpetra_Vector().Update ( 1. , x.castTo<EpetraVectorEngine>().getEpetra_Vector() , 
        1. , y.castTo<EpetraVectorEngine>().getEpetra_Vector() , 0. );
  }

  void EpetraVectorEngine::subtract(const VectorOperations &x, const VectorOperations &y)
  {
    getEpetra_Vector().Update ( 1. , x.castTo<EpetraVectorEngine>().getEpetra_Vector() , 
        -1. , y.castTo<EpetraVectorEngine>().getEpetra_Vector() , 0. );
  }

  void EpetraVectorEngine::multiply(const VectorOperations &x, const VectorOperations &y)
  {
    getEpetra_Vector().Multiply ( 1. , x.castTo<EpetraVectorEngine>().getEpetra_Vector() , 
        y.castTo<EpetraVectorEngine>().getEpetra_Vector() , 0. );
  }

  void EpetraVectorEngine::divide( const VectorOperations &x, const VectorOperations &y)
  {
    getEpetra_Vector().ReciprocalMultiply ( 1. , 
        y.castTo<EpetraVectorEngine>().getEpetra_Vector() , 
        x.castTo<EpetraVectorEngine>().getEpetra_Vector() , 0. );
  }

  void EpetraVectorEngine::reciprocal(const VectorOperations &x)
  {
    getEpetra_Vector().Reciprocal ( x.castTo<EpetraVectorEngine>().getEpetra_Vector() );
  }

  void EpetraVectorEngine::linearSum(double alpha, const VectorOperations &x,
          double beta, const VectorOperations &y)
  {
    getEpetra_Vector().Update ( alpha , x.castTo<EpetraVectorEngine>().getEpetra_Vector() , 
        beta , y.castTo<EpetraVectorEngine>().getEpetra_Vector() , 0. );
  }

  void EpetraVectorEngine::axpy(double alpha, const VectorOperations &x, const VectorOperations &y)
  {
    linearSum ( alpha , x , 1. , y );
  }

  void EpetraVectorEngine::axpby(double alpha, double beta, const VectorOperations &x)
  {
    linearSum ( alpha , x , beta , *this );
  }

  void EpetraVectorEngine::abs(const VectorOperations &x)
  {
    getEpetra_Vector().Abs ( x.castTo<EpetraVectorEngine>().getEpetra_Vector() );
  }

  double EpetraVectorEngine::min(void) const
  {
    double retVal;
    getEpetra_Vector().MinValue ( &retVal );
    return retVal;
  }

  double EpetraVectorEngine::max(void) const
  {
    double retVal;
    getEpetra_Vector().MaxValue ( &retVal );
    return retVal;
  }

  void EpetraVectorEngine::setRandomValues(void)
  {
    getEpetra_Vector().Random();
    abs(*this);
  }


  void EpetraVectorEngine::setValuesByLocalID(int num, int *indices , const double *vals)
  {
    INCREMENT_COUNT("Virtual");
    for ( int i = 0 ; i != num ; i++ )
      getEpetra_Vector()[indices[i]] = vals[i];
  }


  void EpetraVectorEngine::setLocalValuesByGlobalID(int num, int *indices , const double *vals)
  {
    INCREMENT_COUNT("Virtual");
    getEpetra_Vector().ReplaceGlobalValues ( num , const_cast<double *> (vals) , indices );
  }

  void EpetraVectorEngine::addValuesByLocalID(int num, int *indices , const double *vals)
  {
    INCREMENT_COUNT("Virtual");
    for ( int i = 0 ; i != num ; i++ )
      getEpetra_Vector()[indices[i]] += vals[i];
  }

  void EpetraVectorEngine::addLocalValuesByGlobalID(int num, int *indices , const double *vals)
  {
    INCREMENT_COUNT("Virtual");
    getEpetra_Vector().SumIntoGlobalValues ( num , const_cast<double *> (vals) , indices );
  }

  void EpetraVectorEngine::getValuesByLocalID(int num, int *indices , double *vals) const
  {
    INCREMENT_COUNT("Virtual");
    AMP_ERROR( "This shouldn't be called" );
    for ( int i = 0 ; i != num ; i++ )
      vals[i] = getEpetra_Vector()[indices[i]];
  }

  void EpetraVectorEngine::getLocalValuesByGlobalID(int num, int *indices , double *vals) const
  {
    INCREMENT_COUNT("Virtual");
    AMP_ERROR( "This shouldn't be called" );
    for ( int i = 0 ; i != num ; i++ )
      vals[i] = getEpetra_Vector()[indices[i]];
  }
  double EpetraVectorEngine::L1Norm(void) const
  {
    double retVal;
    getEpetra_Vector().Norm1 ( &retVal );
    return retVal;
  }

  double EpetraVectorEngine::L2Norm(void) const
  {
    double retVal;
    getEpetra_Vector().Norm2 ( &retVal );
    return retVal;
  }


  double EpetraVectorEngine::maxNorm(void) const
  {
    double retVal;
    getEpetra_Vector().NormInf ( &retVal );
    return retVal;
  }

  double EpetraVectorEngine::dot(const VectorOperations &x) const
  {
    double retVal;
    getEpetra_Vector().Dot ( x.castTo<EpetraVectorEngine>().getEpetra_Vector() , &retVal );
    return retVal;
  }

}
}

