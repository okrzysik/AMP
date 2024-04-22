#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/CSRMatrix.h"
#include "AMP/matrices/CSRPolicy.h"
#include "AMP/matrices/DenseSerialMatrix.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/data/DenseSerialMatrixData.h"
#include "AMP/utils/Utilities.h"

#ifdef AMP_USE_PETSC
    #include "AMP/matrices/petsc/NativePetscMatrix.h"
    #include "AMP/vectors/petsc/PetscHelpers.h"
#endif
#ifdef AMP_USE_TRILINOS
    #include "AMP/matrices/trilinos/EpetraMatrixData.h"
    #include "AMP/matrices/trilinos/ManagedEpetraMatrix.h"
    #include <Epetra_CrsMatrix.h>
#endif

#include <functional>


namespace AMP::LinearAlgebra {


/********************************************************
 * Check if we have a spare matrix available             *
 ********************************************************/
#if defined( AMP_USE_TRILINOS ) || defined( AMP_USE_PETSC )
bool haveSparseMatrix() { return true; }
#else
bool haveSparseMatrix() { return false; }
#endif


/********************************************************
 * Build a ManagedPetscMatrix                             *
 ********************************************************/
std::shared_ptr<AMP::LinearAlgebra::Matrix>
createManagedMatrix( AMP::LinearAlgebra::Vector::shared_ptr leftVec,
                     AMP::LinearAlgebra::Vector::shared_ptr rightVec,
                     const std::function<std::vector<size_t>( size_t )> &getRow,
                     const std::string &type )
{
    if ( type == "ManagedEpetraMatrix" ) {
#if defined( AMP_USE_TRILINOS )
        // Get the DOFs
        auto leftDOF  = leftVec->getDOFManager();
        auto rightDOF = rightVec->getDOFManager();
        if ( leftDOF->getComm().compare( rightDOF->getComm() ) == 0 )
            AMP_ERROR( "leftDOF and rightDOF on different comm groups is NOT tested" );
        AMP_MPI comm = leftDOF->getComm();
        if ( comm.getSize() == 1 )
            comm = AMP_MPI( AMP_COMM_SELF );

        // Create the matrix parameters
        auto params =
            std::make_shared<AMP::LinearAlgebra::MatrixParameters>( leftDOF, rightDOF, comm );
        params->d_CommListLeft  = leftVec->getCommunicationList();
        params->d_CommListRight = rightVec->getCommunicationList();
        params->d_VariableLeft  = leftVec->getVariable();
        params->d_VariableRight = rightVec->getVariable();

        // Add the row sizes and local columns to the matrix parameters
        std::set<size_t> columns;
        size_t row_start = leftDOF->beginDOF();
        size_t row_end   = leftDOF->endDOF();
        for ( size_t row = row_start; row < row_end; row++ ) {
            auto cols = getRow( row );
            params->setEntriesInRow( row - row_start, cols.size() );
            params->addColumns( cols );
        }

        // Create the matrix
        auto newMatrixData = std::make_shared<AMP::LinearAlgebra::EpetraMatrixData>( params );
        newMatrixData->setEpetraMaps( leftVec, rightVec );
        // Initialize the matrix
        for ( size_t row = row_start; row < row_end; row++ ) {
            auto col = getRow( row );
            newMatrixData->createValuesByGlobalID( row, col );
        }
        newMatrixData->fillComplete();
        auto newMatrix = std::make_shared<AMP::LinearAlgebra::ManagedEpetraMatrix>( newMatrixData );
        newMatrix->zero();
        newMatrix->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
        return newMatrix;
#else
        NULL_USE( leftVec );
        NULL_USE( rightVec );
        NULL_USE( getRow );
        NULL_USE( type );
        AMP_ERROR( "Unable to build ManagedEpetraMatrix without Trilinos" );
#endif
    } else {
        AMP_ERROR( "Unknown ManagedMatrix type" );
        return nullptr;
    }
}

/********************************************************
 * Build a CSRMatrix                                    *
 ********************************************************/
template<typename Policy>
std::shared_ptr<AMP::LinearAlgebra::Matrix>
createCSRMatrix( AMP::LinearAlgebra::Vector::shared_ptr leftVec,
                 AMP::LinearAlgebra::Vector::shared_ptr rightVec,
                 const std::function<std::vector<size_t>( size_t )> &getRow )
{
    // Get the DOFs
    auto leftDOF  = leftVec->getDOFManager();
    auto rightDOF = rightVec->getDOFManager();
    if ( leftDOF->getComm().compare( rightDOF->getComm() ) == 0 )
        AMP_ERROR( "leftDOF and rightDOF on different comm groups is NOT tested, and needs to "
                   "be fixed" );
    AMP_MPI comm = leftDOF->getComm();
    if ( comm.getSize() == 1 )
        comm = AMP_MPI( AMP_COMM_SELF );
    // Create the matrix parameters
    auto params = std::make_shared<AMP::LinearAlgebra::MatrixParameters>( leftDOF, rightDOF, comm );
    params->d_VariableLeft  = leftVec->getVariable();
    params->d_VariableRight = rightVec->getVariable();

    // Add the row sizes and local columns to the matrix parameters
    size_t row_start = leftDOF->beginDOF();
    size_t row_end   = leftDOF->endDOF();
    for ( size_t row = row_start; row < row_end; row++ ) {
        auto cols = getRow( row );
        params->setEntriesInRow( row - row_start, cols.size() );
        params->addColumns( cols );
    }
    params->findUniqueColumns();
    // Create the matrix
    auto data      = std::make_shared<AMP::LinearAlgebra::CSRMatrixData<Policy>>( params );
    auto newMatrix = std::make_shared<AMP::LinearAlgebra::CSRMatrix<Policy>>( data );
    // Initialize the matrix
    newMatrix->zero();
    newMatrix->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
    return newMatrix;
}

using DefaultCSRPolicy = CSRPolicy<size_t, int, double>;

template std::shared_ptr<AMP::LinearAlgebra::Matrix>
createCSRMatrix<DefaultCSRPolicy>( AMP::LinearAlgebra::Vector::shared_ptr leftVec,
                                   AMP::LinearAlgebra::Vector::shared_ptr rightVec,
                                   const std::function<std::vector<size_t>( size_t )> &getRow );


/********************************************************
 * Build a DenseSerialMatrix                             *
 ********************************************************/
std::shared_ptr<AMP::LinearAlgebra::Matrix>
createDenseSerialMatrix( AMP::LinearAlgebra::Vector::shared_ptr leftVec,
                         AMP::LinearAlgebra::Vector::shared_ptr rightVec )
{
    // Get the DOFs
    auto leftDOF  = leftVec->getDOFManager();
    auto rightDOF = rightVec->getDOFManager();
    if ( leftDOF->getComm().compare( rightDOF->getComm() ) == 0 )
        AMP_ERROR( "leftDOF and rightDOF on different comm groups is NOT tested, and needs to "
                   "be fixed" );
    AMP_MPI comm = leftDOF->getComm();
    if ( comm.getSize() == 1 )
        comm = AMP_MPI( AMP_COMM_SELF );
    else
        AMP_ERROR( "serial dense matrix does not support parallel matrices" );
    // Create the matrix parameters
    auto params = std::make_shared<AMP::LinearAlgebra::MatrixParameters>( leftDOF, rightDOF, comm );
    params->d_VariableLeft  = leftVec->getVariable();
    params->d_VariableRight = rightVec->getVariable();
    // Create the matrix
    auto data      = std::make_shared<AMP::LinearAlgebra::DenseSerialMatrixData>( params );
    auto newMatrix = std::make_shared<AMP::LinearAlgebra::DenseSerialMatrix>( data );
    // Initialize the matrix
    newMatrix->zero();
    newMatrix->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
    return newMatrix;
}

/********************************************************
 * Build a NativePetscMatrix                             *
 ********************************************************/
std::shared_ptr<AMP::LinearAlgebra::Matrix>
createNativePetscMatrix( AMP::LinearAlgebra::Vector::shared_ptr leftVec,
                         AMP::LinearAlgebra::Vector::shared_ptr rightVec,
                         const std::function<std::vector<size_t>( size_t )> &getRow )
{
#if defined( AMP_USE_PETSC )
    // Get the DOFs
    auto leftDOF  = leftVec->getDOFManager();
    auto rightDOF = rightVec->getDOFManager();
    if ( leftDOF->getComm().compare( rightDOF->getComm() ) == 0 )
        AMP_ERROR( "leftDOF and rightDOF on different comm groups is NOT tested, and needs to "
                   "be fixed" );
    AMP_MPI comm = leftDOF->getComm();
    // Create the matrix parameters
    auto params = std::make_shared<AMP::LinearAlgebra::MatrixParameters>( leftDOF, rightDOF, comm );
    params->d_VariableLeft  = leftVec->getVariable();
    params->d_VariableRight = rightVec->getVariable();

    // Add the row sizes and local columns to the matrix parameters
    std::set<size_t> columns;
    size_t row_start = leftDOF->beginDOF();
    size_t row_end   = leftDOF->endDOF();
    for ( size_t row = row_start; row < row_end; row++ ) {
        auto cols = getRow( row );
        params->setEntriesInRow( row - row_start, cols.size() );
        params->addColumns( cols );
    }
    // Create the matrix
    auto newMatrix = std::make_shared<AMP::LinearAlgebra::NativePetscMatrix>( params );
    // Initialize the matrix
    //    newMatrix->zero();
    newMatrix->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
    return newMatrix;
#else
    NULL_USE( leftVec );
    NULL_USE( rightVec );
    NULL_USE( getRow );
    AMP_ERROR( "Unable to build NativePetscMatrix without Petsc" );
    return nullptr;
#endif
}


/********************************************************
 * Test the matrix to ensure it is valid                 *
 ********************************************************/
static void test( std::shared_ptr<AMP::LinearAlgebra::Matrix> matrix )
{
    auto leftDOF         = matrix->getLeftDOFManager();
    auto rightDOF        = matrix->getRightDOFManager();
    size_t N_local_row1  = leftDOF->numLocalDOF();
    size_t N_local_row2  = matrix->numLocalRows();
    size_t N_local_col1  = rightDOF->numLocalDOF();
    size_t N_local_col2  = matrix->numLocalColumns();
    size_t N_global_row1 = leftDOF->numGlobalDOF();
    size_t N_global_row2 = matrix->numGlobalRows();
    size_t N_global_col1 = rightDOF->numGlobalDOF();
    size_t N_global_col2 = matrix->numGlobalColumns();
    AMP_ASSERT( N_local_row1 == N_local_row2 );
    AMP_ASSERT( N_local_col1 == N_local_col2 );
    AMP_ASSERT( N_global_row1 == N_global_row2 );
    AMP_ASSERT( N_global_col1 == N_global_col2 );
    AMP_ASSERT( !matrix->getComm().isNull() );
}


/********************************************************
 * Matrix builder                                        *
 ********************************************************/
std::shared_ptr<AMP::LinearAlgebra::Matrix>
createMatrix( AMP::LinearAlgebra::Vector::shared_ptr rightVec,
              AMP::LinearAlgebra::Vector::shared_ptr leftVec,
              const std::string &type,
              std::function<std::vector<size_t>( size_t )> getRow )
{
    // Determine the type of matrix to build
    std::string type2 = type;
    if ( type == "auto" ) {
        // this is only meant to be in the short term
        // once tests can use any matrix we should change
        // exclusively to CSRMatrix
#if defined( AMP_USE_TRILINOS )
        type2 = "ManagedEpetraMatrix";
#else
        type2 = "CSRMatrix";
#endif
    }
    // Create the default getRow function (if not provided)
    if ( !getRow ) {
        const auto leftDOF  = leftVec->getDOFManager().get();
        const auto rightDOF = rightVec->getDOFManager().get();
        getRow              = [leftDOF, rightDOF]( size_t row ) {
            auto elem = leftDOF->getElement( row );
            return rightDOF->getRowDOFs( elem );
        };
    }
    // Build the matrix
    std::shared_ptr<AMP::LinearAlgebra::Matrix> matrix;
    if ( type2 == "ManagedEpetraMatrix" ) {
        matrix = createManagedMatrix( leftVec, rightVec, getRow, type2 );
        test( matrix );
    } else if ( type2 == "NativePetscMatrix" ) {
        matrix = createNativePetscMatrix( leftVec, rightVec, getRow );
        test( matrix );
    } else if ( type2 == "CSRMatrix" ) {
        matrix = createCSRMatrix<DefaultCSRPolicy>( leftVec, rightVec, getRow );
        test( matrix );
    } else if ( type2 == "DenseSerialMatrix" ) {
        matrix = createDenseSerialMatrix( leftVec, rightVec );
        test( matrix );
    } else {
        AMP_ERROR( "Unknown matrix type to build" );
    }
    return matrix;
}


/********************************************************
 * Create Matrix from PETSc Mat                          *
 ********************************************************/
#if defined( AMP_USE_PETSC )
std::shared_ptr<Matrix> createMatrix( Mat M, bool deleteable )
{
    auto matrix = std::make_shared<NativePetscMatrix>( M, deleteable );
    AMP_ASSERT( !matrix->getComm().isNull() );
    return matrix;
}
#endif


} // namespace AMP::LinearAlgebra
