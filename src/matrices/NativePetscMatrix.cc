
//extern "C" {
#include "petscmat.h"
#include "petscvec.h"
//}

#include "Matrix.h"
#include "NativePetscMatrix.h"
#include "vectors/Vector.h"
#include "vectors/NativePetscVector.h"

#include <cassert>


namespace AMP {
namespace LinearAlgebra {

  void NativePetscMatrix::multiply ( shared_ptr other_op , shared_ptr &result )
  {
    if ( !other_op->isA<NativePetscMatrix> () )
    {
      AMP_ERROR( "Incompatible matrix types" );
    }

    NativePetscMatrix *res = new NativePetscMatrix;
    MatMatMult ( d_Mat , other_op->castTo<NativePetscMatrix>().d_Mat , MAT_INITIAL_MATRIX , PETSC_DEFAULT , &(res->d_Mat ));
    result = Matrix::shared_ptr ( res );
  }

  Vector::shared_ptr  NativePetscMatrix::extractDiagonal ( Vector::shared_ptr v )
  {
    Vector::shared_ptr  retVal;
    if ( v->isA<NativePetscVector>() )
    {
      retVal = v;
    }
    else
    {
      retVal = getRightVector();
    }
    MatGetDiagonal ( getMat() , retVal->castTo<PetscVector>().getVec() );
    return retVal;
  }
  
  Vector::shared_ptr NativePetscMatrix::getRightVector ()
  {
    Vec a;
    MatGetVecs ( d_Mat , &a , PETSC_NULL );
    boost::shared_ptr<NativePetscVectorParameters> npvParam ( new NativePetscVectorParameters ( a ) );
    npvParam->d_Deleteable = true;
    return Vector::shared_ptr ( new NativePetscVector ( npvParam ) );
  }

  Vector::shared_ptr NativePetscMatrix::getLeftVector ()
  {
    Vec a;
    MatGetVecs ( d_Mat , PETSC_NULL , &a );
    boost::shared_ptr<NativePetscVectorParameters> npvParam ( new NativePetscVectorParameters ( a ) );
    npvParam->d_Deleteable = true;
    return Vector::shared_ptr ( new NativePetscVector ( npvParam ) );
  }

  void NativePetscMatrix::getRowByGlobalID ( int row , std::vector<unsigned int> &cols , std::vector<double> &values ) const
  {
    int numCols;
    MatGetRow ( d_Mat , row , &numCols , PETSC_NULL , PETSC_NULL );
    cols.resize ( numCols );
    values.resize ( numCols );
    if ( numCols )
    {
      const PetscInt *out_cols;
      const PetscScalar *out_vals;
      MatGetRow ( d_Mat , row , &numCols , &out_cols , &out_vals );
      std::copy ( (unsigned int *)out_cols , (unsigned int *)(out_cols + numCols) , cols.begin() );
      std::copy ( (double *)out_vals , (double *)(out_vals + numCols) , values.begin() );
    }
  }

}
}//end 

