#include "vectors/petsc/NativePetscVector.h"

namespace AMP {
namespace LinearAlgebra {


inline Matrix::shared_ptr  NativePetscMatrix::cloneMatrix () const
{ 
    Mat new_mat;
    MatDuplicate ( d_Mat , MAT_DO_NOT_COPY_VALUES , &new_mat );
    AMP_ERROR( "not qutie implemented" );
    return shared_ptr ( new NativePetscMatrix ( new_mat , true ) );
}


inline size_t NativePetscMatrix::numGlobalRows () const
{
    int rows , cols;
    MatGetSize ( d_Mat , &rows , &cols );
    return rows;
}


inline size_t NativePetscMatrix::numGlobalColumns () const
{
    int rows , cols;
    MatGetSize ( d_Mat , &rows , &cols );
    return cols;
}


inline NativePetscMatrix::NativePetscMatrix ( Mat m , bool internally_created )
{
    d_Mat = m;
    d_MatCreatedInternally = internally_created;
}


inline NativePetscMatrix::NativePetscMatrix() {
    d_MatCreatedInternally = false;
}


inline NativePetscMatrix::~NativePetscMatrix ()
{
    if ( d_MatCreatedInternally ) {
        #if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==0 )
            MatDestroy( d_Mat ); 
        #elif ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==2 )
            MatDestroy( &d_Mat ); 
        #else
            #error Not programmed for this version yet
        #endif
    }
}


inline Matrix::shared_ptr   NativePetscMatrix::duplicateMat ( Mat m )
{
    Mat newMat;
    MatDuplicate ( m , MAT_DO_NOT_COPY_VALUES , &newMat );
    return Matrix::shared_ptr ( new NativePetscMatrix ( newMat , true ) );
}


inline void  NativePetscMatrix::copyFromMat ( Mat m )
{
    MatCopy ( m , d_Mat , SAME_NONZERO_PATTERN );
}


inline void  NativePetscMatrix::mult ( Vector::const_shared_ptr in , Vector::shared_ptr out )
{
    MatMult ( d_Mat , in->castTo<const NativePetscVector>().getVec() , out->castTo<NativePetscVector>().getVec() );
}


inline void  NativePetscMatrix::multTranspose ( Vector::const_shared_ptr in , Vector::shared_ptr out )
{
    MatMultTranspose ( d_Mat , in->castTo<const NativePetscVector>().getVec() , out->castTo<NativePetscVector>().getVec() );
}


inline void  NativePetscMatrix::addValuesByGlobalID ( int num_rows , int num_cols , int *rows, int *cols, double *values )
{
    MatSetValues ( d_Mat , num_rows , rows , num_cols , cols , values , ADD_VALUES );
}


inline void  NativePetscMatrix::setValuesByGlobalID ( int num_rows , int num_cols , int *rows, int *cols, double *values )
{
    MatSetValues ( d_Mat , num_rows , rows , num_cols , cols , values , INSERT_VALUES );
}


inline void  NativePetscMatrix::scale ( double alpha )
{
    MatScale ( d_Mat , alpha );
}


inline void  NativePetscMatrix::axpy ( double alpha , const Matrix  &x )
{
    MatAXPY ( d_Mat , alpha , x.castTo<NativePetscMatrix>().d_Mat , SAME_NONZERO_PATTERN );
}


inline void  NativePetscMatrix::setScalar ( double ans )
{
    if ( ans != 0.0 )
        AMP_ERROR( "Cannot perform operation on NativePetscMatrix yet!" );
    MatZeroEntries ( d_Mat );
}


inline void NativePetscMatrix::setDiagonal ( const Vector::shared_ptr &in )
{
    const PetscVector &pVec = in->castTo<NativePetscVector> ();
    MatDiagonalSet ( d_Mat , pVec.getVec() , INSERT_VALUES );
}
      

inline void NativePetscMatrix::makeConsistent ()
{
    MatAssemblyBegin ( d_Mat , MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd   ( d_Mat , MAT_FINAL_ASSEMBLY );
}


inline double NativePetscMatrix::L1Norm () const
{
    double retVal;
    MatNorm ( d_Mat , NORM_1 , &retVal );
    return retVal;
}


}
}

