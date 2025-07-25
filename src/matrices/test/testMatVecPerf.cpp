#include "AMP/AMP_TPLs.h"
#include "AMP/IO/PIO.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/matrices/CSRConfig.h"
#include "AMP/matrices/CSRMatrix.h"
#include "AMP/matrices/MatrixBuilder.h"
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

#include "ProfilerApp.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

// This test is adapted from testMatVec.cpp and is set up to give some basic
// profiling information regarding matvec products with different matrix
// classes

// Number of products to evaluate to average out timings
#define NUM_PRODUCTS 1000
#define NUM_PRODUCTS_TRANS 100

size_t matVecTestWithDOFs( AMP::UnitTest *ut,
                           std::string type,
                           std::shared_ptr<AMP::Discretization::DOFManager> &dofManager,
                           bool testTranspose,
                           std::string accelerationBackend )
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
#else
    auto inVec     = AMP::LinearAlgebra::createVector( dofManager, inVar );
    auto outVec    = AMP::LinearAlgebra::createVector( dofManager, outVar );
#endif

    // Create the matrix
    auto matrix = AMP::LinearAlgebra::createMatrix(
        inVec, outVec, AMP::Utilities::backendFromString( accelerationBackend ), type );
    if ( matrix ) {
        ut->passes( type + ", " + accelerationBackend + ": Able to create a square matrix" );
    } else {
        ut->failure( type + ", " + accelerationBackend + ": Unable to create a square matrix" );
    }

    fillWithPseudoLaplacian( matrix, dofManager );
    matrix->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );

    size_t nGlobalRows = matrix->numGlobalRows();
    size_t nLocalRows  = matrix->numLocalRows();
    AMP::pout << type << " Global rows: " << nGlobalRows << " Local rows: " << nLocalRows
              << std::endl;

#if defined( AMP_USE_HYPRE )
    using scalar_t = typename AMP::LinearAlgebra::scalar_info<AMP::LinearAlgebra::hypre_real>::type;
#else
    using scalar_t = double;
#endif

    auto x = matrix->getRightVector();
    auto y = matrix->getLeftVector();

    x->setToScalar( 1.0 );
    // this shouldn't be necessary, but evidently is!
    x->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    y->zero();

    for ( int nProd = 0; nProd < NUM_PRODUCTS; ++nProd ) {
        matrix->mult( x, y );
    }

    auto yNorm = static_cast<scalar_t>( y->L1Norm() );

    if ( yNorm == static_cast<scalar_t>( matrix->numGlobalRows() ) ) {
        ut->passes( type + ", " + accelerationBackend +
                    ": Passes 1 norm test with pseudo Laplacian" );
    } else {
        AMP::pout << "1 Norm " << yNorm << ", number of rows " << matrix->numGlobalRows()
                  << std::endl;
        ut->failure( type + ", " + accelerationBackend +
                     ": Fails 1 norm test with pseudo Laplacian" );
    }

    if ( testTranspose && NUM_PRODUCTS_TRANS ) {
        // Repeat test with transpose multiply (Laplacian is symmetric)
        y->setToScalar( 1.0 );
        y->makeConsistent();
        x->zero();
        for ( int nProd = 0; nProd < NUM_PRODUCTS_TRANS; ++nProd ) {
            matrix->multTranspose( y, x );
        }

        auto xNorm = static_cast<scalar_t>( x->L1Norm() );

        if ( xNorm == static_cast<scalar_t>( matrix->numGlobalRows() ) ) {
            ut->passes( type + ", " + accelerationBackend +
                        ": Passes 1 norm test with pseudo Laplacian transpose" );
        } else {
            AMP::pout << "Transpose 1 Norm " << xNorm << ", number of rows "
                      << matrix->numGlobalRows() << std::endl;
            ut->failure( type + ", " + accelerationBackend +
                         ": Fails 1 norm test with pseudo Laplacian transpose" );
        }
    }

    return nGlobalRows;
}

size_t matVecTest( AMP::UnitTest *ut, std::string input_file )
{
    std::string log_file = "output_testMatVecPerf";
    AMP::logOnlyNodeZero( log_file );

    // Read the input file
    auto input_db = AMP::Database::parseInputFile( input_file );

    // Get the Mesh database and create the mesh parameters
    auto database = input_db->getDatabase( "Mesh" );
    auto params   = std::make_shared<AMP::Mesh::MeshParameters>( database );
    auto comm     = AMP::AMP_MPI( AMP_COMM_WORLD );
    params->setComm( comm );

    // Get the acceleration backend for the matrix
    std::vector<std::string> backends;
    if ( input_db->keyExists( "MatrixAccelerationBackend" ) ) {
        backends.emplace_back( input_db->getString( "MatrixAccelerationBackend" ) );
    } else {
        backends.emplace_back( "serial" );
#if ( defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS ) )
        backends.emplace_back( "kokkos" );
#endif
#ifdef USE_DEVICE
        backends.emplace_back( "hip_cuda" );
#endif
    }

    // Create the meshes from the input database
    auto mesh = AMP::Mesh::MeshFactory::create( params );

    // Create the DOF manager
    auto scalarDOFs =
        AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::GeomType::Vertex, 1, 1 );

    // Test on defined matrix types
#if defined( AMP_USE_TRILINOS )
    matVecTestWithDOFs( ut, "ManagedEpetraMatrix", scalarDOFs, true, "serial" );
#endif
#if defined( AMP_USE_PETSC )
    matVecTestWithDOFs( ut, "NativePetscMatrix", scalarDOFs, true, "serial" );
#endif
    size_t nGlobal = 0;
    for ( auto &backend : backends ) {
        nGlobal = matVecTestWithDOFs(
            ut, "CSRMatrix", scalarDOFs, backend == "hip_cuda" ? false : true, backend );
    }
    return nGlobal;
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

        files.emplace_back( "input_testMatVecPerf-1" );
    }

    size_t nGlobal = 0;
    for ( auto &file : files )
        nGlobal = matVecTest( &ut, file );

    ut.report();

    // build unique profile name to avoid collisions
    std::ostringstream ss;
    ss << "testMatVecPerf_r" << std::setw( 3 ) << std::setfill( '0' )
       << AMP::AMPManager::getCommWorld().getSize() << "_n" << std::setw( 9 ) << std::setfill( '0' )
       << nGlobal;
    PROFILE_SAVE( ss.str() );

    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();
    return num_failed;
}
