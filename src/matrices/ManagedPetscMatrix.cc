

#include "petscmat.h"
#include "petscvec.h"

#include "vectors/Vector.h"
#include "vectors/ManagedPetscVector.h"
#include "vectors/trilinos/EpetraVectorEngine.h"

#include "ManagedPetscMatrix.h"



PetscErrorCode  _AMP_Mult ( Mat m , Vec i , Vec o )
{
  void *ctx;
  MatShellGetContext ( m , &ctx );
  AMP::LinearAlgebra::Matrix  *pMatrix = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscMatrix *> ( ctx );
  AMP::LinearAlgebra::Vector  *pVecIn = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *> ( i->data );
  AMP::LinearAlgebra::Vector  *pVecOut = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *> ( o->data );

  AMP::LinearAlgebra::Vector::shared_ptr  pvin ( pVecIn , AMP::LinearAlgebra::ExternalVectorDeleter() );
  AMP::LinearAlgebra::Vector::shared_ptr  pvout ( pVecOut , AMP::LinearAlgebra::ExternalVectorDeleter() );

  pMatrix->mult ( pvin , pvout );
  return 0;
}

PetscErrorCode  _AMP_Mult_add ( Mat m , Vec i , Vec a , Vec o )
{
  void *ctx;
  MatShellGetContext ( m , &ctx );
  AMP::LinearAlgebra::Matrix  *pMatrix = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscMatrix *> ( ctx );
  AMP::LinearAlgebra::Vector  *pVecIn = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *> ( i->data );
  AMP::LinearAlgebra::Vector  *pVecAdd = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *> ( a->data );
  AMP::LinearAlgebra::Vector  *pVecOut = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *> ( o->data );

  AMP::LinearAlgebra::Vector::shared_ptr  pvin ( pVecIn , AMP::LinearAlgebra::ExternalVectorDeleter() );
  AMP::LinearAlgebra::Vector::shared_ptr  pvout ( pVecOut , AMP::LinearAlgebra::ExternalVectorDeleter() );

  pMatrix->mult ( pvin , pvout );
  pVecOut->add ( *pVecOut , *pVecAdd );
  return 0;
}

PetscErrorCode _AMP_GetDiagonal ( Mat m , Vec d )
{
  void *ctx;
  MatShellGetContext ( m , &ctx );

  AMP::LinearAlgebra::Matrix  *pMatrix = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscMatrix *> ( ctx );
  AMP::LinearAlgebra::Vector  *pVecIn = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *> ( d->data );
  AMP::LinearAlgebra::Vector::shared_ptr t ( pVecIn , AMP::LinearAlgebra::ExternalVectorDeleter () );
  pMatrix->extractDiagonal ( t );
  return 0;
}

PetscErrorCode _AMP_GetVecs(Mat m,Vec *right,Vec *left)
{
  void *ctx;
  MatShellGetContext ( m , &ctx );
  AMP::LinearAlgebra::Matrix    *pMatrix = static_cast<AMP::LinearAlgebra::ManagedPetscMatrix *> ( ctx );

  if(right!=PETSC_NULL)
    {
      AMP::LinearAlgebra::Vector::shared_ptr pRight = AMP::LinearAlgebra::PetscVector::view ( pMatrix->getRightVector() );
      VecDuplicate ( pRight->castTo<AMP::LinearAlgebra::PetscVector>().getVec() , right );
    }

  if(left!=PETSC_NULL)
    {
      AMP::LinearAlgebra::Vector::shared_ptr pLeft = AMP::LinearAlgebra::PetscVector::view ( pMatrix->getLeftVector() );
      VecDuplicate ( pLeft->castTo<AMP::LinearAlgebra::PetscVector>().getVec() , left );
    }

  return 0;
}

PetscErrorCode _AMP_AXPY(Mat y,PetscScalar alpha , Mat x , MatStructure )
{
  void *ctx1;
  void *ctx2;
  MatShellGetContext ( x , &ctx1 );
  MatShellGetContext ( y , &ctx2 );
  AMP::LinearAlgebra::Matrix   *pX = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscMatrix *> ( ctx1 );
  AMP::LinearAlgebra::Matrix   *pY = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscMatrix *> ( ctx2 );

  pY->axpy ( alpha , *pX );
  return 0;
}

PetscErrorCode _AMP_Scale (Mat x , PetscScalar alpha )
{
  void *ctx1;
  MatShellGetContext ( x , &ctx1 );
  AMP::LinearAlgebra::Matrix   *pX = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscMatrix *> ( ctx1 );

  pX->scale ( alpha );
  return 0;
}

namespace AMP {
namespace LinearAlgebra {

  void  ManagedPetscMatrix::initPetscMat ()
  {
    ManagedPetscMatrixParameters &params = d_pParameters->castTo<ManagedPetscMatrixParameters> ();
    MPI_Comm petsc_comm = params.getEpetraComm().getCommunicator();
    MatCreateShell ( petsc_comm ,
                     params.getLocalSize(),
                     params.getLocalSize() ,
                     PETSC_DETERMINE ,
                     PETSC_DETERMINE ,
                     static_cast<void *> ( this ) ,
                     &d_Mat );

    MatShellSetOperation ( d_Mat , MATOP_MULT , (void(*)(void)) _AMP_Mult );
    MatShellSetOperation ( d_Mat , MATOP_GET_VECS , (void(*)(void)) _AMP_GetVecs );
    MatShellSetOperation ( d_Mat , MATOP_GET_DIAGONAL , (void(*)(void)) _AMP_GetDiagonal );
    MatShellSetOperation ( d_Mat , MATOP_MULT_ADD , (void(*)(void)) _AMP_Mult_add );
    MatShellSetOperation ( d_Mat , MATOP_AXPY , (void(*)(void)) _AMP_AXPY );
    MatShellSetOperation ( d_Mat , MATOP_SCALE , (void(*)(void)) _AMP_Scale );
  }

  ManagedPetscMatrix::ManagedPetscMatrix ( ParametersPtr params )
    : PetscMatrix ( params )
    , ManagedEpetraMatrix ( params )
  {
//    std::cout << "ManagedPetscMatrix:: WARNING!!!!!! the matrix is currently assumed to be square. This needs to be fixed!!!" << std::endl;
    initPetscMat ();
  }

  ManagedPetscMatrix::ManagedPetscMatrix ( const ManagedPetscMatrix &rhs )
    : Matrix ( rhs.d_pParameters )
    , PetscMatrix ( rhs.d_pParameters )
    , ManagedEpetraMatrix ( rhs )
  {
//    std::cout << "ManagedPetscMatrix:: WARNING!!!!!! the matrix is currently assumed to be square. This needs to be fixed!!!" << std::endl;
    initPetscMat ();
  }

  Matrix::shared_ptr   ManagedPetscMatrix::duplicateMat ( Mat m , AMP_MPI comm )
  {
    int  global_start , global_end , global_size , t;
    MatGetOwnershipRange ( m , &global_start , &global_end );
    MatGetSize ( m , &global_size , &t );
    int  num_local = global_end - global_start;

    ManagedPetscMatrixParameters *params = new ManagedPetscMatrixParameters ( num_local , global_size , global_start , global_size , global_start , comm );
    //ManagedPetscMatrixParameters::iterator  cur_entry = params->begin();
    int i = 0;
    for ( ; global_start != global_end ; global_start++ )
    {
      int num_cols;
      MatGetRow ( m , global_start , &num_cols , PETSC_NULL , PETSC_NULL );
      params->setEntriesInRow ( i , num_cols );
      MatRestoreRow ( m , global_start , &num_cols , PETSC_NULL , PETSC_NULL );
      params->addMapping ( i , global_start );
      i++;
    }
    Matrix::shared_ptr  ret_val ( new ManagedPetscMatrix ( ParametersPtr ( params ) ) );
    return ret_val;
  }

  void  ManagedPetscMatrix::copyFromMat ( Mat m )
  {
    ManagedPetscMatrixParameters::iterator  cur_entry = d_pParameters->castTo<ManagedPetscMatrixParameters>().begin();
    while ( cur_entry != d_pParameters->castTo<ManagedPetscMatrixParameters>().end() )
    {
      int num_cols;
      const int *cols;
      const double  *data;
      MatGetRow ( m , *cur_entry , &num_cols , &cols , &data );
      createValuesByGlobalID ( 1 , num_cols , &*cur_entry , const_cast<int *>(cols) , const_cast<double *>(data) );
      MatRestoreRow ( m , *cur_entry , &num_cols , &cols , &data );
      cur_entry++;
    }
    d_epetraMatrix->FillComplete ();

  }

}
}//end namespace




