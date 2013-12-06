#include "vectors/petsc/NativePetscVector.h"

extern "C"{
#include "assert.h"
}



namespace AMP {
namespace LinearAlgebra {

NativePetscVector::NativePetscVector ( VectorParameters::shared_ptr in_params ): 
    NativeVector(), 
    PetscVector(), 
    VectorEngine()
{ 
    NativePetscVectorParameters &npvParams = in_params->castTo<NativePetscVectorParameters> ();
    d_petscVec = npvParams.d_InVec; 
    d_pArray = 0; 
    CommunicationListParameters::shared_ptr params ( new CommunicationListParameters () );
    params->d_comm = npvParams.d_Comm;
    params->d_localsize = npvParams.d_localsize;
    setCommunicationList ( CommunicationList::shared_ptr ( new CommunicationList ( params ) ) );
    d_bDeleteMe = npvParams.d_Deleteable;
    d_DOFManager = AMP::Discretization::DOFManager::shared_ptr( new AMP::Discretization::DOFManager( npvParams.d_localsize, npvParams.d_Comm ) );
}


NativePetscVector::~NativePetscVector ()
{
    resetArray();
    if ( d_bDeleteMe ) {
        #if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==0 )
            VecDestroy(d_petscVec);
        #elif ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==2 )
            VecDestroy(&d_petscVec);
        #else
            #error Not programmed for this version yet
        #endif
    }
}


Vector::shared_ptr NativePetscVector::cloneVector(const Variable::shared_ptr var ) const 
{ 
    resetArray();
    Vec  new_petscVec;
    VecDuplicate ( d_petscVec , &new_petscVec );
    boost::shared_ptr<NativePetscVectorParameters> npvParams( new NativePetscVectorParameters( new_petscVec, true ) );
    npvParams->d_Comm = getComm();
    Vector::shared_ptr retVal =  Vector::shared_ptr ( new NativePetscVector ( npvParams ) );
    retVal->setVariable ( var );
    return retVal;
}

void NativePetscVector::putRawData ( const double *in )
{
    int a , b;
    VecGetOwnershipRange ( d_petscVec , &a , &b );
    int *offs = new int [b-a];
    for ( int j = 0 ; j != b-a ; j++ )
      offs[j] = a+j;
    VecSetValues ( d_petscVec , b-a , offs , in , INSERT_VALUES );
}


void NativePetscVector::copyOutRawData ( double *out ) const
{
    std::copy ( getRawDataBlock<double> ( 0 ) ,
                getRawDataBlock<double> ( 0 ) + getLocalSize() ,
                out );
}


}
}

