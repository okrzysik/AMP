#ifdef USE_AMP_VECTORS

#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/DenseSerialMatrix.h"
#include "AMP/utils/Utilities.h"

#ifdef USE_EXT_TRILINOS
#include "AMP/matrices/trilinos/ManagedEpetraMatrix.h"
#include "AMP/vectors/trilinos/epetra/EpetraVectorEngine.h"
#ifdef USE_EXT_PETSC
#include "AMP/matrices/petsc/ManagedPetscMatrix.h"
#include "AMP/vectors/petsc/ManagedPetscVector.h"
#endif
#endif

#include <functional>


namespace AMP {
namespace LinearAlgebra {


/********************************************************
 * Build a ManagedPetscMatrix                             *
 ********************************************************/
AMP::LinearAlgebra::Matrix::shared_ptr
createManagedMatrix( AMP::LinearAlgebra::Vector::shared_ptr leftVec,
                     AMP::LinearAlgebra::Vector::shared_ptr rightVec,
                     const std::function<std::vector<size_t>( size_t )> &getRow,
                     const std::string &type )
{
#if defined( USE_EXT_TRILINOS )
    // Get the DOFs
    auto leftDOF  = leftVec->getDOFManager();
    auto rightDOF = rightVec->getDOFManager();
    if ( leftDOF->getComm().compare( rightVec->getComm() ) == 0 )
        AMP_ERROR( "leftDOF and rightDOF on different comm groups is NOT tested, and needs to "
                   "be fixed" );
    AMP_MPI comm = leftDOF->getComm();
    if ( comm.getSize() == 1 )
        comm = AMP_MPI( AMP_COMM_SELF );

    // Create the matrix parameters
    auto params = std::make_shared<AMP::LinearAlgebra::ManagedEpetraMatrixParameters>(
        leftDOF, rightDOF, comm );
    params->d_CommListLeft  = leftVec->getCommunicationList();
    params->d_CommListRight = rightVec->getCommunicationList();
    params->d_VariableLeft  = leftVec->getVariable();
    params->d_VariableRight = rightVec->getVariable();

    // Add the row sizes and local columns to the matrix parameters
    std::set<size_t> columns;
    size_t row_start = leftDOF->beginDOF();
    size_t row_end   = leftDOF->endDOF();
    for ( size_t row = row_start; row < row_end; row++ ) {
        auto col = getRow( row );
        params->setEntriesInRow( row - row_start, col.size() );
        for ( auto &tmp : col )
            columns.insert( tmp );
    }
    params->addColumns( columns );

    // Create the matrix
    std::shared_ptr<AMP::LinearAlgebra::ManagedEpetraMatrix> newMatrix;
    if ( type == "ManagedPetscMatrix" ) {
#if defined( USE_EXT_PETSC )
        newMatrix.reset( new AMP::LinearAlgebra::ManagedPetscMatrix( params ) );
#else
        AMP_ERROR( "Unable to build ManagedPetscMatrix without PETSc" );
#endif
    } else if ( type == "ManagedEpetraMatrix" ) {
        newMatrix.reset( new AMP::LinearAlgebra::ManagedEpetraMatrix( params ) );
    } else {
        AMP_ERROR( "Unknown ManagedMatrix type" );
    }

    // Initialize the matrix
    for ( size_t row = row_start; row < row_end; row++ ) {
        auto col = getRow( row );
        newMatrix->createValuesByGlobalID( row, col );
    }
    std::dynamic_pointer_cast<AMP::LinearAlgebra::EpetraMatrix>( newMatrix )
        ->setEpetraMaps( leftVec, rightVec );
    newMatrix->fillComplete();
    newMatrix->zero();
    newMatrix->makeConsistent();
    return newMatrix;
#else
    NULL_USE( leftVec );
    NULL_USE( rightVec );
    NULL_USE( type );
    NULL_USE( getRow );
    AMP_ERROR( "Unable to build a ManagedMatrix without TRILINOS" );
    return AMP::LinearAlgebra::Matrix::shared_ptr();
#endif
}


/********************************************************
 * Build a DenseSerialMatrix                             *
 ********************************************************/
AMP::LinearAlgebra::Matrix::shared_ptr
createDenseSerialMatrix( AMP::LinearAlgebra::Vector::shared_ptr leftVec,
                         AMP::LinearAlgebra::Vector::shared_ptr rightVec )
{
    // Get the DOFs
    auto leftDOF  = leftVec->getDOFManager();
    auto rightDOF = rightVec->getDOFManager();
    if ( leftDOF->getComm().compare( rightVec->getComm() ) == 0 )
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
    auto newMatrix = std::make_shared<AMP::LinearAlgebra::DenseSerialMatrix>( params );
    // Initialize the matrix
    newMatrix->zero();
    newMatrix->makeConsistent();
    return newMatrix;
}


/********************************************************
 * Test the matrix to ensure it is valid                 *
 ********************************************************/
static void test( AMP::LinearAlgebra::Matrix::shared_ptr matrix )
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
}


/********************************************************
 * Matrix builder                                        *
 ********************************************************/
AMP::LinearAlgebra::Matrix::shared_ptr
createMatrix( AMP::LinearAlgebra::Vector::shared_ptr rightVec,
              AMP::LinearAlgebra::Vector::shared_ptr leftVec,
              const std::string &type,
              std::function<std::vector<size_t>( size_t )> getRow )
{
    // Determine the type of matrix to build
    std::string type2 = type;
    if ( type == "auto" ) {
#if defined( USE_EXT_TRILINOS ) && defined( USE_EXT_PETSC )
        type2 = "ManagedPetscMatrix";
#elif defined( USE_EXT_TRILINOS )
        type2 = "ManagedEpetraMatrix";
#else
        type2 = "DenseSerialMatrix";
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
    AMP::LinearAlgebra::Matrix::shared_ptr matrix;
    if ( type2 == "ManagedPetscMatrix" || type2 == "ManagedEpetraMatrix" ) {
        matrix = createManagedMatrix( leftVec, rightVec, getRow, type2 );
    } else if ( type2 == "DenseSerialMatrix" ) {
        matrix = createDenseSerialMatrix( leftVec, rightVec );
    } else {
        AMP_ERROR( "Unknown matrix type to build" );
    }
    // Run some quick checks on the matrix
    if ( matrix )
        test( matrix );
    return matrix;
}
} // namespace LinearAlgebra
} // namespace AMP

#endif
