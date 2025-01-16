#include "AMP/AMP_TPLs.h"
#include "AMP/IO/PIO.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/matrices/CSRMatrix.h"
#include "AMP/matrices/CSRMatrixParameters.h"
#include "AMP/matrices/CSRPolicy.h"
#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/testHelpers/MatrixDataTransforms.h"
#include "AMP/matrices/testHelpers/MatrixTests.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#if defined( AMP_USE_HYPRE )
    #include "AMP/matrices/data/hypre/HypreCSRPolicy.h"
#endif

#include "ProfilerApp.h"

#include <iomanip>
#include <memory>
#include <sstream>
#include <string>

// This test is adapted from testMatVecPerf.cpp (and hence testMatVec.cpp)
// This version forms simple aggregates of nodes of a matrix (A) into another
// matrix (P) then tests that A*(P*x)=(A*P)*x for some random vector x

// TODO: extend this to triple product once transposing is supported

template<typename Policy, typename Allocator>
void findAggregates( std::shared_ptr<AMP::LinearAlgebra::CSRLocalMatrixData<Policy, Allocator>> A,
                     std::vector<typename Policy::lidx_t> &agg_id,
                     typename Policy::lidx_t &num_agg )
{
    using lidx_t = typename Policy::lidx_t;

    const auto A_nrows = static_cast<lidx_t>( A->numLocalRows() );

    // vectors for aggregate associations and aggregate sizes
    agg_id.resize( A_nrows, -1 ); // -1 for not assn
    std::vector<lidx_t> agg_size( A_nrows, 0 );

    // unpack members of A_diag
    auto [A_rs, A_cols, A_cols_loc, A_coeffs] = A->getDataFields();

    // first pass initilizes aggregates from nodes that have no
    // neighbors that are already associated
    num_agg = 0;
    for ( lidx_t row = 0; row < A_nrows; ++row ) {
        // check if any members of this row are already associated
        bool have_assn = false;
        for ( lidx_t c = A_rs[row]; c < A_rs[row + 1]; ++c ) {
            have_assn = have_assn || ( agg_id[A_cols_loc[c]] >= 0 );
        }
        // if no associations create new aggregate from row
        if ( !have_assn ) {
            for ( lidx_t c = A_rs[row]; c < A_rs[row + 1]; ++c ) {
                agg_id[A_cols_loc[c]] = num_agg;
                agg_size[num_agg]++;
            }
            // increment current id to start working on next aggregate
            ++num_agg;
        }
    }

    // second pass adds unmarked entries to the smallest aggregate they are nbrs with
    // entries are unmarked because they neigbored some aggregate in the above,
    // thus every unmarked entry will neighbor some aggregate
    for ( lidx_t row = 0; row < A_nrows; ++row ) {
        if ( agg_id[row] >= 0 ) {
            // this row already assigned, skip ahead
            continue;
        }
        // find smallest neighboring aggregate
        lidx_t small_agg_id = -1, small_agg_size = A_nrows + 1;
        for ( lidx_t c = A_rs[row]; c < A_rs[row + 1]; ++c ) {
            const auto id = agg_id[A_cols_loc[c]];
            if ( id >= 0 && ( agg_size[id] < small_agg_size ) ) {
                small_agg_size = agg_size[id];
                small_agg_id   = id;
            }
        }
        AMP_ASSERT( small_agg_id >= 0 );
        agg_id[row] = small_agg_id;
        agg_size[small_agg_id]++;
    }
}

template<typename Policy, typename Allocator>
std::shared_ptr<AMP::LinearAlgebra::CSRMatrix<Policy, Allocator>>
createAggregateMatrix( std::shared_ptr<AMP::LinearAlgebra::CSRMatrix<Policy, Allocator>> A )
{
    using lidx_t = typename Policy::lidx_t;

    // get aggregates just using diagonal block
    std::vector<lidx_t> agg_id;
    lidx_t num_agg;
    auto A_data = std::dynamic_pointer_cast<AMP::LinearAlgebra::CSRMatrixData<Policy, Allocator>>(
        A->getMatrixData() );
    findAggregates( A_data->getDiagMatrix(), agg_id, num_agg );

    // Build matrix parameters object
    auto leftDOFs      = A->getRightDOFManager(); // inner dof manager for A*P
    auto rightDOFs     = std::make_shared<AMP::Discretization::DOFManager>( num_agg, A->getComm() );
    auto leftClParams  = std::make_shared<AMP::LinearAlgebra::CommunicationListParameters>();
    auto rightClParams = std::make_shared<AMP::LinearAlgebra::CommunicationListParameters>();
    leftClParams->d_comm         = A->getComm();
    rightClParams->d_comm        = A->getComm();
    leftClParams->d_localsize    = leftDOFs->numLocalDOF();
    rightClParams->d_localsize   = rightDOFs->numLocalDOF();
    leftClParams->d_remote_DOFs  = leftDOFs->getRemoteDOFs();
    rightClParams->d_remote_DOFs = rightDOFs->getRemoteDOFs();

    // Create parameters, variables, and data
    auto params = std::make_shared<AMP::LinearAlgebra::MatrixParameters>(
        leftDOFs, rightDOFs, A->getComm(), A_data->getRightVariable(), A_data->getRightVariable() );
    auto P_data = std::make_shared<AMP::LinearAlgebra::CSRMatrixData<Policy, Allocator>>( params );

    // non-zeros only in diag block and only one per row
    P_data->setNNZ( std::vector<lidx_t>( agg_id.size(), 1 ),
                    std::vector<lidx_t>( agg_id.size(), 0 ) );

    // fill in data (diag block only) using aggregates from above
    auto P_diag                               = P_data->getDiagMatrix();
    auto [P_rs, P_cols, P_cols_loc, P_coeffs] = P_diag->getDataFields();
    const auto begin_col                      = rightDOFs->beginDOF();
    const auto num_row                        = static_cast<lidx_t>( leftDOFs->numLocalDOF() );
    for ( lidx_t row = 0; row < num_row; ++row ) {
        const auto agg = agg_id[row];
        const auto rs  = P_rs[row];
        P_cols[rs]     = begin_col + agg;
        P_cols_loc[rs] = agg;
        P_coeffs[rs]   = 1.0;
    }

    // reset dof managers and return matrix
    P_data->resetDOFManagers();
    return std::make_shared<AMP::LinearAlgebra::CSRMatrix<Policy, Allocator>>( P_data );
}

size_t matMultTestWithDOFs( AMP::UnitTest *ut,
                            std::shared_ptr<AMP::Discretization::DOFManager> &dofManager )
{
    auto comm = AMP::AMP_MPI( AMP_COMM_WORLD );
    // Create the vectors
    auto inVar  = std::make_shared<AMP::LinearAlgebra::Variable>( "inputVar" );
    auto outVar = std::make_shared<AMP::LinearAlgebra::Variable>( "outputVar" );
#ifdef USE_DEVICE
    auto inVec = AMP::LinearAlgebra::createVector(
        dofManager, inVar, true, AMP::Utilities::MemoryType::managed );
    auto outVec = AMP::LinearAlgebra::createVector(
        dofManager, outVar, true, AMP::Utilities::MemoryType::managed );
    ut->expected_failure( "SpGEMM for device CSRMatrix not implemented" );
    return 0;
#else
    auto inVec             = AMP::LinearAlgebra::createVector( dofManager, inVar );
    auto outVec            = AMP::LinearAlgebra::createVector( dofManager, outVar );
#endif

    // Create the matrix
    auto A = AMP::LinearAlgebra::createMatrix( inVec, outVec, "CSRMatrix" );
    fillWithPseudoLaplacian( A, dofManager );
    A->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );

    size_t nGlobalRows = A->numGlobalRows();
    size_t nLocalRows  = A->numLocalRows();
    AMP::pout << "CSRMatrix Global rows: " << nGlobalRows << " Local rows: " << nLocalRows
              << std::endl;

#if defined( AMP_USE_HYPRE )
    using DefaultCSRPolicy = AMP::LinearAlgebra::CSRPolicy<HYPRE_BigInt, HYPRE_Int, HYPRE_Real>;
#else
    using DefaultCSRPolicy = AMP::LinearAlgebra::CSRPolicy<size_t, int, double>;
#endif

    // Create aggregate matrix
    auto Acsr = std::dynamic_pointer_cast<
        AMP::LinearAlgebra::CSRMatrix<DefaultCSRPolicy, AMP::HostAllocator<void>>>( A );
    auto P = createAggregateMatrix( Acsr );

    // perform A*P SpGEMM
    auto AP = AMP::LinearAlgebra::Matrix::matMultiply( A, P );

    // vectors of ones to apply operators to
    auto xa  = A->getRightVector();
    auto xp  = P->getRightVector();
    auto xap = AP->getRightVector();
    xa->setToScalar( 1.0 );
    xp->setToScalar( 1.0 );
    xap->setToScalar( 1.0 );
    xa->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    xp->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    xap->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );

    // output vectors
    auto ya  = A->getLeftVector();
    auto yp  = AP->getLeftVector();
    auto yap = AP->getLeftVector();
    ya->zero();
    yp->zero();
    yap->zero();

    // vector of ones is in col-span of P so A * xa should be same as AP * xap
    A->mult( xa, ya );
    P->mult( xp, yp );
    AP->mult( xap, yap );
    const auto l1ya    = static_cast<double>( ya->L1Norm() );
    const auto l1yp    = static_cast<double>( yp->L1Norm() );
    const auto l1yap   = static_cast<double>( yap->L1Norm() );
    const auto nrows_d = static_cast<double>( nGlobalRows );

    if ( AMP::Utilities::approx_equal( l1ya, l1yap ) &&
         AMP::Utilities::approx_equal( l1yp, nrows_d ) ) {
        ut->passes( "matMultiply A*P CSRMatrix" );
    } else {
        AMP::pout << "matMultiply A*P CSRMatrix fails with l1ya = " << l1ya << ", l1yp = " << l1yp
                  << ", l1yap = " << l1yap << std::endl;
        ut->failure( "matMultiply A*P CSRMatrix" );
    }

    return nGlobalRows;
}

size_t matMultTest( AMP::UnitTest *ut, std::string input_file )
{
    std::string log_file = "output_testMatMultCoarsen";
    AMP::logOnlyNodeZero( log_file );

    // Read the input file
    auto input_db = AMP::Database::parseInputFile( input_file );

    // Get the Mesh database and create the mesh parameters
    auto database = input_db->getDatabase( "Mesh" );
    auto params   = std::make_shared<AMP::Mesh::MeshParameters>( database );
    auto comm     = AMP::AMP_MPI( AMP_COMM_WORLD );
    params->setComm( comm );

    // Create the meshes from the input database
    auto mesh = AMP::Mesh::MeshFactory::create( params );

    // Create the DOF manager
    auto scalarDOFs =
        AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::GeomType::Vertex, 1, 1 );

    return matMultTestWithDOFs( ut, scalarDOFs );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    std::vector<std::string> files;
    PROFILE_ENABLE();

    if ( argc > 1 ) {

        files.emplace_back( argv[1] );

    } else {
        // Same input files as testMatVecPerf
        files.emplace_back( "input_testMatVecPerf-1" );
    }

    size_t nGlobal = 0;
    for ( auto &file : files )
        nGlobal = matMultTest( &ut, file );

    ut.report();

    // build unique profile name to avoid collisions
    std::ostringstream ss;
    ss << "testMatMultCoarsen_r" << std::setw( 3 ) << std::setfill( '0' )
       << AMP::AMPManager::getCommWorld().getSize() << "_n" << std::setw( 9 ) << std::setfill( '0' )
       << nGlobal;
    PROFILE_SAVE( ss.str() );

    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();
    return num_failed;
}
