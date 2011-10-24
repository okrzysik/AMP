#include "NativePetscVector.h"

extern "C"{
#include "assert.h"
}



namespace AMP {
namespace LinearAlgebra {

  NativePetscVector::NativePetscVector ( VectorParameters::shared_ptr in_params ) 
    : NativeVector() , PetscVector () , VectorEngine ()
  { 
    NativePetscVectorParameters &npvParams = in_params->castTo<NativePetscVectorParameters> ();
    d_petscVec = npvParams.d_InVec; 
    d_pArray = 0; 
    CommunicationListParameters::shared_ptr params ( new CommunicationListParameters () );
    params->d_comm = npvParams.d_Comm;
    setCommunicationList ( CommunicationList::shared_ptr ( new CommunicationList ( params ) ) );
    d_bDeleteMe = npvParams.d_Deleteable;
  }

  NativePetscVector::~NativePetscVector ()
  {
    resetArray();
    if ( d_bDeleteMe )
    {
      VecDestroy ( d_petscVec );
    }
  }

  Vector::shared_ptr NativePetscVector::cloneVector(const Variable::shared_ptr var ) const 
  { 
    resetArray();
    Vec  new_petscVec;
    VecDuplicate ( d_petscVec , &new_petscVec );
    boost::shared_ptr<NativePetscVectorParameters> npvParams ( new NativePetscVectorParameters ( new_petscVec ) );
    npvParams->d_Comm = getComm();
    npvParams->d_Deleteable = true;
    Vector::shared_ptr retVal =  Vector::shared_ptr ( new NativePetscVector ( npvParams ) );
    retVal->setVariable ( var );
    return retVal;
  }

  void NativePetscVector::putRawData ( double *in )
  {
    int a , b;
    VecGetOwnershipRange ( d_petscVec , &a , &b );
    int *offs = new int [b-a];
    for ( int j = 0 ; j != b-a ; j++ )
      offs[j] = a+j;
    VecSetValues ( d_petscVec , b-a , offs , in , INSERT_VALUES );
  }


}
}

