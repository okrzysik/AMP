#include "utils/Utilities.h"

#include "EpetraVectorEngine.h"

#ifdef USE_MPI
    #include <Epetra_MpiComm.h>
#else
    #include <Epetra_SerialComm.h>
#endif

namespace AMP {
namespace LinearAlgebra {


/********************************************************
* EpetraVectorEngineParameters constructors             *
********************************************************/
EpetraVectorEngineParameters::EpetraVectorEngineParameters ( size_t local_size , size_t global_size , AMP_MPI c ):
    VectorEngineParameters(),
    d_comm( c.getCommunicator() ) 
{
    d_global = global_size;
    d_comm.sumScan( &local_size, &d_end, 1 );
    d_begin = d_end - local_size;
}
EpetraVectorEngineParameters::EpetraVectorEngineParameters ( size_t local_size , size_t global_size , boost::shared_ptr<Epetra_Map> emap , AMP_MPI ecomm ):
    VectorEngineParameters(),
    d_emap( emap ),
    d_comm( ecomm )
{
    d_global = global_size;
    d_comm.sumScan( &local_size, &d_end, 1 );
    d_begin = d_end - local_size;
}


/********************************************************
* Function to return (and create) the Epetra_Map        *
********************************************************/
Epetra_Map  &EpetraVectorEngineParameters::getEpetraMap()
{
    if ( d_emap.get() != NULL )
        return *d_emap;
    // Create the epetra map
    #ifdef USE_MPI
        Epetra_MpiComm  comm = d_comm.getCommunicator();
    #else
        Epetra_SerialComm  comm;
    #endif
    AMP_INSIST(d_global<0x80000000,"Epetra does not support vectors with global size greater than 2^31");
    size_t local_size = d_end-d_begin;
    //std::vector<int> ids(local_size,0);
    //for (size_t i=0; i<local_size; i++)
    //    ids[i] = (int) (i+d_begin);
    //d_emap = boost::shared_ptr<Epetra_Map> ( new Epetra_Map ( (int) d_global, (int) local_size, &ids[0], 0, comm ) );
    d_emap = boost::shared_ptr<Epetra_Map> ( new Epetra_Map ( (int) d_global, (int) local_size, 0, comm ) );
    // Check the map to make sure it is correct
    AMP_ASSERT(local_size==(size_t)d_emap->NumMyPoints());
    AMP_ASSERT(d_global==(size_t)d_emap->NumGlobalPoints());
    AMP_ASSERT(d_begin==(size_t)d_emap->MinMyGID());
    AMP_ASSERT(d_end-1==(size_t)d_emap->MaxMyGID());
    AMP_ASSERT(0==(size_t)d_emap->MinAllGID());
    AMP_ASSERT(d_global-1==(size_t)d_emap->MaxAllGID());
    AMP_ASSERT(0==(size_t)d_emap->MinLID());
    AMP_ASSERT(local_size-1==(size_t)d_emap->MaxLID());
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


  void EpetraVectorEngine::setValuesByLocalID(int num, size_t *indices , const double *vals)
  {
    INCREMENT_COUNT("Virtual");
    for ( int i = 0 ; i != num ; i++ )
      getEpetra_Vector()[indices[i]] = vals[i];
  }


  void EpetraVectorEngine::setLocalValuesByGlobalID(int num, size_t *indices , const double *vals)
  {
    INCREMENT_COUNT("Virtual");
    if ( num==0 )
        return;
    AMP_ASSERT(getGlobalSize()<0x80000000);
    std::vector<int> indices2(num,0);
    for (int i=0; i<num; i++)
        indices2[i] = (int) indices[i];
    getEpetra_Vector().ReplaceGlobalValues ( num , const_cast<double *> (vals) , &indices2[0] );
  }

  void EpetraVectorEngine::addValuesByLocalID(int num, size_t *indices , const double *vals)
  {
    INCREMENT_COUNT("Virtual");
    if ( num==0 )
        return;
    for ( int i = 0 ; i != num ; i++ )
      getEpetra_Vector()[indices[i]] += vals[i];
  }

  void EpetraVectorEngine::addLocalValuesByGlobalID(int num, size_t *indices , const double *vals)
  {
    INCREMENT_COUNT("Virtual");
    if ( num==0 )
        return;
    AMP_ASSERT(getGlobalSize()<0x80000000);
    std::vector<int> indices2(num,0);
    for (int i=0; i<num; i++)
        indices2[i] = (int) indices[i];
    getEpetra_Vector().SumIntoGlobalValues ( num , const_cast<double *> (vals) , &indices2[0] );
  }

  void EpetraVectorEngine::getValuesByLocalID(int num, size_t *indices , double *vals) const
  {
    INCREMENT_COUNT("Virtual");
    AMP_ERROR( "This shouldn't be called" );
    //for ( int i = 0 ; i != num ; i++ )
    //  vals[i] = getEpetra_Vector()[indices[i]];
  }

  void EpetraVectorEngine::getLocalValuesByGlobalID(int num, size_t *indices , double *vals) const
  {
    INCREMENT_COUNT("Virtual");
    AMP_ERROR( "This shouldn't be called" );
    //for ( int i = 0 ; i != num ; i++ )
    //  vals[i] = getEpetra_Vector()[indices[i]];
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

