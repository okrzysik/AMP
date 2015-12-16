
// extern "C" {
#include "petscmat.h"
#include "petscvec.h"
//}

#include "matrices/Matrix.h"
#include "matrices/petsc/NativePetscMatrix.h"
#include "vectors/Vector.h"
#include "vectors/petsc/NativePetscVector.h"

#include <cassert>


namespace AMP {
namespace LinearAlgebra {


/********************************************************
* Multiply the matrix by a vector                       *
********************************************************/
void NativePetscMatrix::multiply( shared_ptr other_op, shared_ptr &result )
{
    if ( !other_op->isA<NativePetscMatrix>() ) {
        AMP_ERROR( "Incompatible matrix types" );
    }

    auto res = new NativePetscMatrix;
    MatMatMult( d_Mat,
                other_op->castTo<NativePetscMatrix>().d_Mat,
                MAT_INITIAL_MATRIX,
                PETSC_DEFAULT,
                &( res->d_Mat ) );
    result = Matrix::shared_ptr( res );
}


/********************************************************
* Get the diagonal                                      *
********************************************************/
Vector::shared_ptr NativePetscMatrix::extractDiagonal( Vector::shared_ptr v ) const
{
    Vector::shared_ptr retVal;
    if ( v->isA<NativePetscVector>() ) {
        retVal = v;
    } else {
        retVal = getRightVector();
    }
    MatGetDiagonal( getMat(), retVal->castTo<PetscVector>().getVec() );
    return retVal;
}


/********************************************************
* Get the left/right Vector/DOFManager                  *
********************************************************/
Vector::shared_ptr NativePetscMatrix::getRightVector() const
{
    Vec a;
    MatGetVecs( d_Mat, &a, PETSC_NULL );
    AMP::shared_ptr<NativePetscVectorParameters> npvParam(
        new NativePetscVectorParameters( a, true ) );
    return Vector::shared_ptr( new NativePetscVector( npvParam ) );
}
Vector::shared_ptr NativePetscMatrix::getLeftVector() const
{
    Vec a;
    MatGetVecs( d_Mat, PETSC_NULL, &a );
    AMP::shared_ptr<NativePetscVectorParameters> npvParam(
        new NativePetscVectorParameters( a, true ) );
    return Vector::shared_ptr( new NativePetscVector( npvParam ) );
}
Discretization::DOFManager::shared_ptr NativePetscMatrix::getRightDOFManager() const
{
    return getRightVector()->getDOFManager();
}
Discretization::DOFManager::shared_ptr NativePetscMatrix::getLeftDOFManager() const
{
    return getLeftVector()->getDOFManager();
}


/********************************************************
* Get values/row by global id                           *
********************************************************/
void NativePetscMatrix::getValuesByGlobalID(
    int num_rows, int num_cols, int *rows, int *cols, double *values ) const
{
    // Zero out the data in values
    for ( int i   = 0; i < num_rows * num_cols; i++ )
        values[i] = 0.0;
    // Get the data for each row
    Discretization::DOFManager::shared_ptr leftDOFManager = getLeftDOFManager();
    int firstRow                                          = leftDOFManager->beginDOF();
    int numRows                                           = leftDOFManager->endDOF();
    for ( int i = 0; i < num_rows; i++ ) {
        if ( rows[i] < firstRow || rows[i] >= firstRow + numRows )
            continue;
        int numCols = 0;
        MatGetRow( d_Mat, rows[i], &numCols, PETSC_NULL, PETSC_NULL );
        if ( numCols == 0 )
            continue;
        const PetscInt *out_cols;
        const PetscScalar *out_vals;
        MatGetRow( d_Mat, rows[i], &numCols, &out_cols, &out_vals );
        for ( int j1 = 0; j1 < num_cols; j1++ ) {
            for ( int j2 = 0; j2 < numCols; j2++ ) {
                if ( cols[j1] == out_cols[j2] )
                    values[i * num_cols + j1] = out_vals[j2];
            }
        }
    }
}
void NativePetscMatrix::getRowByGlobalID( int row,
                                          std::vector<unsigned int> &cols,
                                          std::vector<double> &values ) const
{
    int numCols;
    MatGetRow( d_Mat, row, &numCols, PETSC_NULL, PETSC_NULL );
    cols.resize( numCols );
    values.resize( numCols );
    if ( numCols ) {
        const PetscInt *out_cols;
        const PetscScalar *out_vals;
        MatGetRow( d_Mat, row, &numCols, &out_cols, &out_vals );
        std::copy(
            (unsigned int *) out_cols, (unsigned int *) ( out_cols + numCols ), cols.begin() );
        std::copy( (double *) out_vals, (double *) ( out_vals + numCols ), values.begin() );
    }
}
}
} // end
