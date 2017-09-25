#include "vectors/petsc/NativePetscVector.h"

namespace AMP {
namespace LinearAlgebra {


inline Matrix::shared_ptr NativePetscMatrix::cloneMatrix() const
{
    Mat new_mat;
    MatDuplicate( d_Mat, MAT_DO_NOT_COPY_VALUES, &new_mat );
    AMP_ERROR( "not quite implemented" );
    return shared_ptr( new NativePetscMatrix( new_mat, true ) );
}


inline size_t NativePetscMatrix::numGlobalRows() const
{
    int rows, cols;
    MatGetSize( d_Mat, &rows, &cols );
    return (size_t) rows;
}


inline size_t NativePetscMatrix::numGlobalColumns() const
{
    int rows, cols;
    MatGetSize( d_Mat, &rows, &cols );
    return (size_t) cols;
}


inline NativePetscMatrix::NativePetscMatrix( Mat m, bool internally_created )
{
    d_Mat                  = m;
    d_MatCreatedInternally = internally_created;
}


inline NativePetscMatrix::NativePetscMatrix() { d_MatCreatedInternally = false; }


inline NativePetscMatrix::~NativePetscMatrix()
{
    if ( d_MatCreatedInternally )
        PETSC::matDestroy( &d_Mat );
}


inline Matrix::shared_ptr NativePetscMatrix::duplicateMat( Mat m )
{
    Mat newMat;
    MatDuplicate( m, MAT_DO_NOT_COPY_VALUES, &newMat );
    return Matrix::shared_ptr( new NativePetscMatrix( newMat, true ) );
}


inline void NativePetscMatrix::copyFromMat( Mat m ) { MatCopy( m, d_Mat, SAME_NONZERO_PATTERN ); }


inline void NativePetscMatrix::mult( Vector::const_shared_ptr in, Vector::shared_ptr out )
{
    MatMult( d_Mat,
             dynamic_pointer_cast<const NativePetscVector>( in )->getVec(),
             dynamic_pointer_cast<NativePetscVector>( out )->getVec() );
}


inline void NativePetscMatrix::multTranspose( Vector::const_shared_ptr in, Vector::shared_ptr out )
{
    MatMultTranspose( d_Mat,
                      dynamic_pointer_cast<const NativePetscVector>( in )->getVec(),
                      dynamic_pointer_cast<NativePetscVector>( out )->getVec() );
}


inline void NativePetscMatrix::addValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values )
{
    std::vector<PetscInt> petsc_rows( num_rows );
    std::vector<PetscInt> petsc_cols( num_cols );
    std::copy( rows, rows + num_rows, petsc_rows.begin() );
    std::copy( cols, cols + num_cols, petsc_cols.begin() );

    MatSetValues( d_Mat, num_rows, &petsc_rows[0], num_cols, &petsc_cols[0], values, ADD_VALUES );
}


inline void NativePetscMatrix::setValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values )
{
    std::vector<PetscInt> petsc_rows( num_rows );
    std::vector<PetscInt> petsc_cols( num_cols );
    std::copy( rows, rows + num_rows, petsc_rows.begin() );
    std::copy( cols, cols + num_cols, petsc_cols.begin() );

    MatSetValues(
        d_Mat, num_rows, &petsc_rows[0], num_cols, &petsc_cols[0], values, INSERT_VALUES );
}


inline void NativePetscMatrix::scale( double alpha ) { MatScale( d_Mat, alpha ); }


inline void NativePetscMatrix::axpy( double alpha, const Matrix &x )
{
    AMP_ASSERT( x.numGlobalRows() == this->numGlobalRows() );
    AMP_ASSERT( x.numGlobalColumns() == this->numGlobalColumns() );
    MatAXPY(
        d_Mat, alpha, dynamic_cast<const NativePetscMatrix *>( &x )->d_Mat, SAME_NONZERO_PATTERN );
}


inline void NativePetscMatrix::setScalar( double ans )
{
    if ( ans != 0.0 )
        AMP_ERROR( "Cannot perform operation on NativePetscMatrix yet!" );
    MatZeroEntries( d_Mat );
}


inline void NativePetscMatrix::zero() { MatZeroEntries( d_Mat ); }


inline void NativePetscMatrix::setDiagonal( Vector::const_shared_ptr in )
{
    auto pVec = dynamic_pointer_cast<const NativePetscVector>( in );
    MatDiagonalSet( d_Mat, pVec->getVec(), INSERT_VALUES );
}

inline void NativePetscMatrix::setIdentity()
{
    MatZeroEntries( d_Mat );
    MatShift( d_Mat, 1.0 );
}


inline void NativePetscMatrix::makeConsistent()
{
    MatAssemblyBegin( d_Mat, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( d_Mat, MAT_FINAL_ASSEMBLY );
}


inline double NativePetscMatrix::L1Norm() const
{
    double retVal;
    MatNorm( d_Mat, NORM_1, &retVal );
    return retVal;
}
} // namespace LinearAlgebra
} // namespace AMP
