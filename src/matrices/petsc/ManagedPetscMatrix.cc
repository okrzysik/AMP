#include "petscmat.h"
#include "petscvec.h"

#include "AMP/vectors/ExternalVectorDeleter.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/petsc/ManagedPetscVector.h"
#include "AMP/vectors/petsc/PetscVector.h"
#include "AMP/vectors/trilinos/epetra/EpetraVectorEngine.h"

#include "AMP/matrices/petsc/ManagedPetscMatrix.h"
#include "AMP/matrices/trilinos/ManagedEpetraMatrix.h"


PetscErrorCode _AMP_Mult( Mat m, Vec i, Vec o )
{
    void *ctx;
    MatShellGetContext( m, &ctx );
    AMP::LinearAlgebra::Matrix *pMatrix =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscMatrix *>( ctx );
    AMP::LinearAlgebra::Vector *pVecIn =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( i->data );
    AMP::LinearAlgebra::Vector *pVecOut =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( o->data );

    AMP::LinearAlgebra::Vector::shared_ptr pvin( pVecIn,
                                                 AMP::LinearAlgebra::ExternalVectorDeleter() );
    AMP::LinearAlgebra::Vector::shared_ptr pvout( pVecOut,
                                                  AMP::LinearAlgebra::ExternalVectorDeleter() );

    pMatrix->mult( pvin, pvout );
    return 0;
}


PetscErrorCode _AMP_Mult_add( Mat m, Vec i, Vec a, Vec o )
{
    void *ctx;
    MatShellGetContext( m, &ctx );
    AMP::LinearAlgebra::Matrix *pMatrix =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscMatrix *>( ctx );
    AMP::LinearAlgebra::Vector *pVecIn =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( i->data );
    AMP::LinearAlgebra::Vector *pVecAdd =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( a->data );
    AMP::LinearAlgebra::Vector *pVecOut =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( o->data );

    AMP::LinearAlgebra::Vector::shared_ptr pvin( pVecIn,
                                                 AMP::LinearAlgebra::ExternalVectorDeleter() );
    AMP::LinearAlgebra::Vector::shared_ptr pvout( pVecOut,
                                                  AMP::LinearAlgebra::ExternalVectorDeleter() );

    pMatrix->mult( pvin, pvout );
    pVecOut->add( *pVecOut, *pVecAdd );
    return 0;
}


PetscErrorCode _AMP_GetDiagonal( Mat m, Vec d )
{
    void *ctx;
    MatShellGetContext( m, &ctx );

    AMP::LinearAlgebra::Matrix *pMatrix =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscMatrix *>( ctx );
    AMP::LinearAlgebra::Vector *pVecIn =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( d->data );
    AMP::LinearAlgebra::Vector::shared_ptr t( pVecIn, AMP::LinearAlgebra::ExternalVectorDeleter() );
    pMatrix->extractDiagonal( t );
    return 0;
}


PetscErrorCode _AMP_GetVecs( Mat m, Vec *right, Vec *left )
{
    void *ctx;
    MatShellGetContext( m, &ctx );
    AMP::LinearAlgebra::Matrix *pMatrix =
        static_cast<AMP::LinearAlgebra::ManagedPetscMatrix *>( ctx );
    if ( right != PETSC_NULL ) {
        AMP::shared_ptr<AMP::LinearAlgebra::PetscVector> pRight =
            AMP::dynamic_pointer_cast<AMP::LinearAlgebra::PetscVector>(
                AMP::LinearAlgebra::PetscVector::view( pMatrix->getRightVector() ) );
        VecDuplicate( pRight->getVec(), right );
    }
    if ( left != PETSC_NULL ) {
        AMP::shared_ptr<AMP::LinearAlgebra::PetscVector> pLeft =
            AMP::dynamic_pointer_cast<AMP::LinearAlgebra::PetscVector>(
                AMP::LinearAlgebra::PetscVector::view( pMatrix->getLeftVector() ) );
        VecDuplicate( pLeft->getVec(), left );
    }
    return 0;
}


PetscErrorCode _AMP_AXPY( Mat y, PetscScalar alpha, Mat x, MatStructure )
{
    void *ctx1;
    void *ctx2;
    MatShellGetContext( x, &ctx1 );
    MatShellGetContext( y, &ctx2 );
    AMP::LinearAlgebra::Matrix *pX =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscMatrix *>( ctx1 );
    AMP::LinearAlgebra::Matrix *pY =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscMatrix *>( ctx2 );

    pY->axpy( alpha, *pX );
    return 0;
}


PetscErrorCode _AMP_Scale( Mat x, PetscScalar alpha )
{
    void *ctx1;
    MatShellGetContext( x, &ctx1 );
    AMP::LinearAlgebra::Matrix *pX =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscMatrix *>( ctx1 );

    pX->scale( alpha );
    return 0;
}


namespace AMP {
namespace LinearAlgebra {

using ManagedPetscMatrixParameters = ManagedEpetraMatrixParameters;


void ManagedPetscMatrix::initPetscMat()
{
    AMP::shared_ptr<ManagedPetscMatrixParameters> params =
        AMP::dynamic_pointer_cast<ManagedPetscMatrixParameters>( d_pParameters );
    MPI_Comm petsc_comm = params->getEpetraComm().getCommunicator();
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

    MatShellSetOperation( d_Mat, MATOP_MULT, (void ( * )()) _AMP_Mult );
    MatShellSetOperation( d_Mat, MATOP_GET_VECS, (void ( * )()) _AMP_GetVecs );
    MatShellSetOperation( d_Mat, MATOP_GET_DIAGONAL, (void ( * )()) _AMP_GetDiagonal );
    MatShellSetOperation( d_Mat, MATOP_MULT_ADD, (void ( * )()) _AMP_Mult_add );
    MatShellSetOperation( d_Mat, MATOP_AXPY, (void ( * )()) _AMP_AXPY );
    MatShellSetOperation( d_Mat, MATOP_SCALE, (void ( * )()) _AMP_Scale );
}


ManagedPetscMatrix::ManagedPetscMatrix( MatrixParameters::shared_ptr params )
    : Matrix( params ),
      PetscMatrix( params ),
      ManagedEpetraMatrix( AMP::dynamic_pointer_cast<ManagedEpetraMatrixParameters>( params ) )
{
    //    std::cout << "ManagedPetscMatrix:: WARNING!!!!!! the matrix is currently assumed to be
    //    square. This needs to
    //    be fixed!!!" << std::endl;
    initPetscMat();
}


ManagedPetscMatrix::ManagedPetscMatrix( const ManagedPetscMatrix &rhs )
    : Matrix( rhs.d_pParameters ), PetscMatrix( rhs.d_pParameters ), ManagedEpetraMatrix( rhs )
{
    //    std::cout << "ManagedPetscMatrix:: WARNING!!!!!! the matrix is currently assumed to be
    //    square. This needs to
    //    be fixed!!!" << std::endl;
    initPetscMat();
}


Matrix::shared_ptr ManagedPetscMatrix::duplicateMat( Mat m, AMP_MPI comm )
{
    int global_start, global_end, global_size, t;
    MatGetOwnershipRange( m, &global_start, &global_end );
    MatGetSize( m, &global_size, &t );
    int num_local = global_end - global_start;

    AMP::Discretization::DOFManager::shared_ptr left(
        new AMP::Discretization::DOFManager( num_local, comm ) );
    AMP::Discretization::DOFManager::shared_ptr right(
        new AMP::Discretization::DOFManager( num_local, comm ) );
    ManagedPetscMatrixParameters *params = new ManagedPetscMatrixParameters( left, right, comm );
    // ManagedPetscMatrixParameters::iterator  cur_entry = params->begin();
    int i = 0;
    for ( ; global_start != global_end; global_start++ ) {
        int num_cols;
        MatGetRow( m, global_start, &num_cols, PETSC_NULL, PETSC_NULL );
        params->setEntriesInRow( i, num_cols );
        MatRestoreRow( m, global_start, &num_cols, PETSC_NULL, PETSC_NULL );
        i++;
    }
    Matrix::shared_ptr ret_val( new ManagedPetscMatrix( MatrixParameters::shared_ptr( params ) ) );
    return ret_val;
}


void ManagedPetscMatrix::copyFromMat( Mat m )
{
    AMP::shared_ptr<ManagedPetscMatrixParameters> params =
        AMP::dynamic_pointer_cast<ManagedPetscMatrixParameters>( d_pParameters );
    AMP::Discretization::DOFManager::shared_ptr rowDOF = params->getLeftDOFManager();
    // AMP::Discretization::DOFManager::shared_ptr colDOF = params->getRightDOFManager();
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
    d_epetraMatrix->FillComplete();
}


ManagedPetscMatrix::~ManagedPetscMatrix() { PETSC::matDestroy( &d_Mat ); }


Matrix::shared_ptr ManagedPetscMatrix::cloneMatrix() const
{
    return shared_ptr( new ManagedPetscMatrix( *this ) );
}
} // namespace LinearAlgebra
} // namespace AMP
