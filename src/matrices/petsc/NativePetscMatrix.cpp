#include "petscmat.h"
#include "petscvec.h"

#include "AMP/matrices/Matrix.h"
#include "AMP/matrices/petsc/NativePetscMatrix.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/petsc/NativePetscVector.h"


namespace AMP {
namespace LinearAlgebra {


/********************************************************
 * Multiply the matrix by a vector                       *
 ********************************************************/
void NativePetscMatrix::multiply( shared_ptr other_op, shared_ptr &result )
{
    auto other = std::dynamic_pointer_cast<NativePetscMatrix>( other_op );
    AMP_INSIST( other != nullptr, "Incompatible matrix types" );

    auto res = new NativePetscMatrix;
    MatMatMult( d_Mat, other->d_Mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &( res->d_Mat ) );
    result = Matrix::shared_ptr( res );
}


/********************************************************
 * Get the diagonal                                      *
 ********************************************************/
Vector::shared_ptr NativePetscMatrix::extractDiagonal( Vector::shared_ptr v ) const
{
    Vector::shared_ptr retVal;
    if ( std::dynamic_pointer_cast<NativePetscVector>( v ) ) {
        retVal = v;
    } else {
        retVal = getRightVector();
    }
    MatGetDiagonal( getMat(), std::dynamic_pointer_cast<PetscVector>( retVal )->getVec() );
    return retVal;
}


/********************************************************
 * Get the left/right Vector/DOFManager                  *
 ********************************************************/
Vector::shared_ptr NativePetscMatrix::getRightVector() const
{
    Vec a;
#if PETSC_VERSION_LE( 3, 2, 0 )
    MatGetVecs( d_Mat, &a, PETSC_NULL );
#else
    MatCreateVecs( d_Mat, &a, PETSC_NULL );
#endif
    return std::make_shared<NativePetscVector>( a, true );
}
Vector::shared_ptr NativePetscMatrix::getLeftVector() const
{
    Vec a;
#if PETSC_VERSION_LE( 3, 2, 0 )
    MatGetVecs( d_Mat, PETSC_NULL, &a );
#else
    MatCreateVecs( d_Mat, PETSC_NULL, &a );
#endif
    return std::make_shared<NativePetscVector>( a, true );
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
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) const
{
    // Zero out the data in values
    for ( size_t i = 0; i < num_rows * num_cols; i++ )
        values[i] = 0.0;
    // Get the data for each row
    Discretization::DOFManager::shared_ptr leftDOFManager = getLeftDOFManager();
    size_t firstRow                                       = leftDOFManager->beginDOF();
    size_t numRows                                        = leftDOFManager->endDOF();

    for ( size_t i = 0; i < num_rows; i++ ) {
        if ( rows[i] < firstRow || rows[i] >= firstRow + numRows )
            continue;
        int numCols = 0;
        MatGetRow( d_Mat, rows[i], &numCols, PETSC_NULL, PETSC_NULL );
        if ( numCols == 0 )
            continue;
        const PetscInt *out_cols;
        const PetscScalar *out_vals;
        MatGetRow( d_Mat, rows[i], &numCols, &out_cols, &out_vals );
        for ( size_t j1 = 0; j1 < num_cols; j1++ ) {
            for ( int j2 = 0; j2 < numCols; j2++ ) {
                if ( (int) cols[j1] == out_cols[j2] )
                    values[i * num_cols + j1] = out_vals[j2];
            }
        }
    }
}
void NativePetscMatrix::getRowByGlobalID( size_t row,
                                          std::vector<size_t> &cols,
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

std::vector<size_t> NativePetscMatrix::getColumnIDs( size_t row ) const
{
    int numCols;
    MatGetRow( d_Mat, row, &numCols, PETSC_NULL, PETSC_NULL );
    std::vector<size_t> cols( numCols );

    if ( numCols ) {
        const PetscInt *out_cols;
        const PetscScalar *out_vals;
        MatGetRow( d_Mat, row, &numCols, &out_cols, &out_vals );
        std::copy(
            (unsigned int *) out_cols, (unsigned int *) ( out_cols + numCols ), cols.begin() );
    }

    return cols;
}
} // namespace LinearAlgebra
} // namespace AMP
