#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/CSRMatrix.h"
#include "AMP/matrices/CSRPolicy.h"
#include "AMP/matrices/DenseSerialMatrix.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/data/DenseSerialMatrixData.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/memory.h"

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
 * Build a ManagedEpetraMatrix                           *
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
        auto params = std::make_shared<AMP::LinearAlgebra::MatrixParameters>(
            leftDOF, rightDOF, comm, getRow );
        params->d_CommListLeft  = leftVec->getCommunicationList();
        params->d_CommListRight = rightVec->getCommunicationList();
        params->d_VariableLeft  = leftVec->getVariable();
        params->d_VariableRight = rightVec->getVariable();

        // Create the matrix
        auto newMatrixData = std::make_shared<AMP::LinearAlgebra::EpetraMatrixData>( params );
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
template<typename Policy, class Allocator>
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
    auto params =
        std::make_shared<AMP::LinearAlgebra::MatrixParameters>( leftDOF, rightDOF, comm, getRow );
    params->d_CommListLeft  = leftVec->getCommunicationList();
    params->d_CommListRight = rightVec->getCommunicationList();
    params->d_VariableLeft  = leftVec->getVariable();
    params->d_VariableRight = rightVec->getVariable();

    // Create the matrix
    auto data = std::make_shared<AMP::LinearAlgebra::CSRMatrixData<Policy, Allocator>>( params );
    auto newMatrix = std::make_shared<AMP::LinearAlgebra::CSRMatrix<Policy, Allocator>>( data );
    // Initialize the matrix
    newMatrix->zero();
    newMatrix->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
    return newMatrix;
}

using DefaultCSRPolicy = CSRPolicy<size_t, int, double>;

template std::shared_ptr<AMP::LinearAlgebra::Matrix>
createCSRMatrix<DefaultCSRPolicy, AMP::HostAllocator<int>>(
    AMP::LinearAlgebra::Vector::shared_ptr leftVec,
    AMP::LinearAlgebra::Vector::shared_ptr rightVec,
    const std::function<std::vector<size_t>( size_t )> &getRow );

#ifdef USE_DEVICE
template std::shared_ptr<AMP::LinearAlgebra::Matrix>
createCSRMatrix<DefaultCSRPolicy, AMP::ManagedAllocator<int>>(
    AMP::LinearAlgebra::Vector::shared_ptr leftVec,
    AMP::LinearAlgebra::Vector::shared_ptr rightVec,
    const std::function<std::vector<size_t>( size_t )> &getRow );

template std::shared_ptr<AMP::LinearAlgebra::Matrix>
createCSRMatrix<DefaultCSRPolicy, AMP::DeviceAllocator<int>>(
    AMP::LinearAlgebra::Vector::shared_ptr leftVec,
    AMP::LinearAlgebra::Vector::shared_ptr rightVec,
    const std::function<std::vector<size_t>( size_t )> &getRow );
#endif


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
    auto params =
        std::make_shared<AMP::LinearAlgebra::MatrixParameters>( leftDOF, rightDOF, comm, getRow );
    params->d_VariableLeft  = leftVec->getVariable();
    params->d_VariableRight = rightVec->getVariable();

    // Create the matrix
    auto newMatrix = std::make_shared<AMP::LinearAlgebra::NativePetscMatrix>( params );
    // Initialize the matrix
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
              std::string type,
              std::function<std::vector<size_t>( size_t )> getRow )
{
    if ( type == "auto" )
        type = DEFAULT_MATRIX; // Definition set by CMake variable DEFAULT_MATRIX
    if ( type == "NULL" )
        return nullptr; // Special case to return nullptr
    // Create the default getRow function (if not provided)
    if ( !getRow ) {
        const auto leftDOF  = leftVec->getDOFManager().get();
        const auto rightDOF = rightVec->getDOFManager().get();
        getRow              = [leftDOF, rightDOF]( size_t row ) {
            auto elem = leftDOF->getElement( row );
            return rightDOF->getRowDOFs( elem );
        };
    }

    // Find memory type associated with (right) vector
    auto memType = AMP::Utilities::getMemoryType( rightVec->getRawDataBlockAsVoid( 0 ) );

    // Build the matrix
    std::shared_ptr<AMP::LinearAlgebra::Matrix> matrix;
    if ( type == "ManagedEpetraMatrix" ) {
        matrix = createManagedMatrix( leftVec, rightVec, getRow, type );
    } else if ( type == "NativePetscMatrix" ) {
        matrix = createNativePetscMatrix( leftVec, rightVec, getRow );
    } else if ( type == "CSRMatrix" ) {
        if ( memType <= AMP::Utilities::MemoryType::host ) {
            AMP::pout << "Building csrmatrix in host memory" << std::endl;
            matrix = createCSRMatrix<DefaultCSRPolicy, AMP::HostAllocator<int>>(
                leftVec, rightVec, getRow );
        } else if ( memType == AMP::Utilities::MemoryType::managed ) {
#ifdef USE_DEVICE
            AMP::pout << "Building csrmatrix in managed memory" << std::endl;
            matrix = createCSRMatrix<DefaultCSRPolicy, AMP::ManagedAllocator<int>>(
                leftVec, rightVec, getRow );
#else
            AMP_ERROR( "Creating CSRMatrix in managed memory requires HIP or CUDA support" );
#endif
        } else if ( memType == AMP::Utilities::MemoryType::device ) {
#ifdef USE_DEVICE
            AMP::pout << "Building csrmatrix in device memory" << std::endl;
            matrix = createCSRMatrix<DefaultCSRPolicy, AMP::DeviceAllocator<int>>(
                leftVec, rightVec, getRow );
#else
            AMP_ERROR( "Creating CSRMatrix in device memory requires HIP or CUDA support" );
#endif
        } else {
            AMP_ERROR( "Unknown memory space in createMatrix" );
        }
    } else if ( type == "DenseSerialMatrix" ) {
        matrix = createDenseSerialMatrix( leftVec, rightVec );
    } else {
        AMP_ERROR( "Unknown matrix type to build" );
    }
    // Check that the matrix is valid
    test( matrix );
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
