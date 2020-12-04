#ifdef USE_EXT_TRILINOS

#include "petsc/private/vecimpl.h"
#include "petscmat.h"
#include "petscvec.h"

#include "AMP/matrices/ManagedMatrixParameters.h"
#include "AMP/matrices/petsc/ManagedPetscMatrix.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/petsc/PetscHelpers.h"
#include "AMP/vectors/petsc/PetscVector.h"


#if PETSC_VERSION_LT( 3, 7, 5 )
#error AMP only supports PETSc 3.7.5 or greater
#endif


PetscErrorCode _AMP_Mult( Mat m, Vec i, Vec o )
{
    auto mat    = PETSC::getAMP( m );
    auto vecIn  = PETSC::getAMP( i );
    auto vecOut = PETSC::getAMP( o );
    mat->mult( vecIn, vecOut );
    return 0;
}


PetscErrorCode _AMP_Mult_add( Mat m, Vec i, Vec a, Vec o )
{
    auto mat    = PETSC::getAMP( m );
    auto vecIn  = PETSC::getAMP( i );
    auto vecAdd = PETSC::getAMP( a );
    auto vecOut = PETSC::getAMP( o );
    mat->mult( vecIn, vecOut );
    vecOut->add( *vecOut, *vecAdd );
    return 0;
}


PetscErrorCode _AMP_GetDiagonal( Mat m, Vec d )
{
    auto mat = PETSC::getAMP( m );
    auto vec = PETSC::getAMP( d );
    mat->extractDiagonal( vec );
    return 0;
}


PetscErrorCode _AMP_GetVecs( Mat m, Vec *right, Vec *left )
{
    auto mat = PETSC::getAMP( m );
    if ( right != PETSC_NULL ) {
        auto pRight = AMP::LinearAlgebra::PetscVector::view( mat->getRightVector() );
        VecDuplicate( pRight->getVec(), right );
    }
    if ( left != PETSC_NULL ) {
        auto pLeft = AMP::LinearAlgebra::PetscVector::view( mat->getLeftVector() );
        VecDuplicate( pLeft->getVec(), left );
    }
    return 0;
}


PetscErrorCode _AMP_AXPY( Mat y, PetscScalar alpha, Mat x, MatStructure )
{
    auto matX = PETSC::getAMP( x );
    auto matY = PETSC::getAMP( y );
    matY->axpy( alpha, *matX );
    return 0;
}


PetscErrorCode _AMP_Scale( Mat x, PetscScalar alpha )
{
    auto mat = PETSC::getAMP( x );
    mat->scale( alpha );
    return 0;
}


namespace AMP {
namespace LinearAlgebra {


void ManagedPetscMatrix::initPetscMat()
{
    auto params         = std::dynamic_pointer_cast<ManagedMatrixParameters>( d_pParameters );
    MPI_Comm petsc_comm = params->getComm().getCommunicator();
    size_t N_col_local  = params->getLocalNumberOfColumns();
    size_t N_row_local  = params->getLocalNumberOfRows();
    size_t N_col_global = params->getGlobalNumberOfColumns();
    size_t N_row_global = params->getGlobalNumberOfRows();

    MatCreateShell( petsc_comm,
                    N_row_local,
                    N_col_local,
                    N_row_global,
                    N_col_global,
                    static_cast<void *>( this ),
                    &d_Mat );

#if PETSC_VERSION_GE( 3, 12, 2 )
    // BP: have not checked for versions above 3.7.5
    MatShellSetManageScalingShifts( d_Mat );
#endif
    MatShellSetOperation( d_Mat, MATOP_MULT, (void ( * )()) _AMP_Mult );
#if PETSC_VERSION_GE( 3, 12, 2 )
    MatShellSetOperation( d_Mat, MATOP_CREATE_VECS, (void ( * )()) _AMP_GetVecs );
#elif PETSC_VERSION_LE( 3, 7, 5 )
    MatShellSetOperation( d_Mat, MATOP_GET_VECS, (void ( * )()) _AMP_GetVecs );
#else
#error Not programmed for this version of petsc
#endif
    MatShellSetOperation( d_Mat, MATOP_GET_DIAGONAL, (void ( * )()) _AMP_GetDiagonal );
    MatShellSetOperation( d_Mat, MATOP_MULT_ADD, (void ( * )()) _AMP_Mult_add );
    MatShellSetOperation( d_Mat, MATOP_AXPY, (void ( * )()) _AMP_AXPY );
    MatShellSetOperation( d_Mat, MATOP_SCALE, (void ( * )()) _AMP_Scale );
}


ManagedPetscMatrix::ManagedPetscMatrix( MatrixParameters::shared_ptr params )
    : Matrix( params ),
      PetscMatrix( params ),
      ManagedEpetraMatrix( std::dynamic_pointer_cast<ManagedMatrixParameters>( params ) )
{
    initPetscMat();
}


ManagedPetscMatrix::ManagedPetscMatrix( const ManagedPetscMatrix &rhs )
    : Matrix( rhs.d_pParameters ), PetscMatrix( rhs.d_pParameters ), ManagedEpetraMatrix( rhs )
{
    initPetscMat();
}


Matrix::shared_ptr ManagedPetscMatrix::duplicateMat( Mat m, AMP_MPI comm )
{
    int global_start, global_end, global_size, t;
    MatGetOwnershipRange( m, &global_start, &global_end );
    MatGetSize( m, &global_size, &t );
    int num_local = global_end - global_start;

    auto left   = std::make_shared<AMP::Discretization::DOFManager>( num_local, comm );
    auto right  = std::make_shared<AMP::Discretization::DOFManager>( num_local, comm );
    auto params = std::make_shared<ManagedMatrixParameters>( left, right, comm );
    int i       = 0;
    for ( ; global_start != global_end; global_start++ ) {
        int num_cols;
        MatGetRow( m, global_start, &num_cols, PETSC_NULL, PETSC_NULL );
        params->setEntriesInRow( i, num_cols );
        MatRestoreRow( m, global_start, &num_cols, PETSC_NULL, PETSC_NULL );
        i++;
    }
    auto ret_val = std::make_shared<ManagedPetscMatrix>( params );
    return ret_val;
}


void ManagedPetscMatrix::copyFromMat( Mat m )
{
    auto params = std::dynamic_pointer_cast<ManagedMatrixParameters>( d_pParameters );
    auto rowDOF = params->getLeftDOFManager();
    for ( size_t i = rowDOF->beginDOF(); i < rowDOF->endDOF(); i++ ) {
        AMP_ASSERT( i < 0x80000000 ); // We have not converted matrices to 64-bits yet
        auto row = (int) i;
        int num_cols;
        const int *cols;
        const double *data;
        MatGetRow( m, row, &num_cols, &cols, &data );
        std::vector<size_t> cols2( num_cols );
        for ( int i = 0; i < num_cols; i++ )
            cols2[i] = cols[i];
        createValuesByGlobalID( row, cols2 );
        MatRestoreRow( m, row, &num_cols, &cols, &data );
    }
    FillComplete();
}


ManagedPetscMatrix::~ManagedPetscMatrix() { PETSC::matDestroy( &d_Mat ); }


Matrix::shared_ptr ManagedPetscMatrix::cloneMatrix() const
{
    return std::shared_ptr<ManagedPetscMatrix>( new ManagedPetscMatrix( *this ) );
}
} // namespace LinearAlgebra
} // namespace AMP

#endif
