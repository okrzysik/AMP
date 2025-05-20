#include "AMP/AMP_TPLs.h"
#include "AMP/IO/PIO.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/solvers/SolverFactory.h"
#include "AMP/solvers/testHelpers/SolverTestParameters.h"
#include "AMP/solvers/testHelpers/testSolverHelpers.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/VectorBuilder.h"

#include <iomanip>
#include <memory>
#include <string>

void linearThermalTest( AMP::UnitTest *ut,
                        const std::string &inputFileName,
                        std::string &accelerationBackend )
{
    PROFILE( "DRIVER::linearThermalTest" );

    // Input and output file names
    std::string input_file = inputFileName;
    std::ostringstream ss;
    ss << "output_testLinSolveRobin_r" << std::setw( 3 ) << std::setfill( '0' )
       << AMP::AMPManager::getCommWorld().getSize();

    AMP::pout << "Running linearThermalTest with input " << input_file << " with "
              << accelerationBackend << " backend" << std::endl;

    // Fill the database from the input file.
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Print from all cores into the output files
    AMP::logAllNodes( ss.str() );

    // Create the Mesh
    const auto meshAdapter = createMesh( input_db );

    auto PowerInWattsVec = constructNeutronicsPowerSource( input_db, meshAdapter );

    // Set appropriate acceleration backend
    auto op_db = input_db->getDatabase( "DiffusionBVPOperator" );
    op_db->putScalar( "AccelerationBackend", accelerationBackend );
#ifdef USE_DEVICE
    op_db->putScalar( "MemoryLocation", "managed" );
#endif
    // Create the Thermal BVP Operator
    auto diffusionOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "DiffusionBVPOperator", input_db ) );

    auto linearOp               = diffusionOperator->getVolumeOperator();
    auto TemperatureInKelvinVec = linearOp->getLeftVector();
    auto RightHandSideVec       = linearOp->getRightVector();

    auto boundaryOpCorrectionVec = RightHandSideVec->clone();

    // Add the boundary conditions corrections
    auto boundaryOp = diffusionOperator->getBoundaryOperator();
    boundaryOp->addRHScorrection( boundaryOpCorrectionVec );
    RightHandSideVec->subtract( *PowerInWattsVec, *boundaryOpCorrectionVec );

    auto &comm = meshAdapter->getComm();

    auto linearSolver = AMP::Solver::Test::buildSolver(
        "LinearSolver", input_db, comm, nullptr, diffusionOperator );

    // Set initial guess
    TemperatureInKelvinVec->setToScalar( 1.0 );

    AMP::pout << "System size: " << RightHandSideVec->getGlobalSize() << std::endl;

    // Use a random initial guess?
    linearSolver->setZeroInitialGuess( false );

    // Solve the problem.
    {
        PROFILE( "DRIVER::linearThermalTest(solve call)" );
        linearSolver->apply( RightHandSideVec, TemperatureInKelvinVec );
    }

    checkConvergence( linearSolver.get(), input_db, input_file, *ut );
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

        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-CG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-IPCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-FCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-CylMesh-CG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-GMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-FGMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BiCGSTAB" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-TFQMR" );

        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-CG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-IPCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-FCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-CG-FCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-GMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-FGMRES" );
        files.emplace_back(
            "input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-GMRESR-GMRES" );
        files.emplace_back(
            "input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-GMRESR-BiCGSTAB" );
        files.emplace_back(
            "input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-GMRESR-TFQMR" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-BiCGSTAB" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-TFQMR" );
#ifdef AMP_USE_PETSC
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-PetscCG" );
        files.emplace_back(
            "input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-PetscFGMRES" );
        //        files.emplace_back(
        //        "input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-PetscBiCGSTAB" );
#endif

#ifdef AMP_USE_HYPRE
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-IPCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG-FCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG-CG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FGMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-GCR" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-GMRES" );
        files.emplace_back(
            "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-BiCGSTAB" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-TFQMR" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-BiCGSTAB" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-TFQMR" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-HypreCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreGMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreGMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-HypreGMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreBiCGSTAB" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreBiCGSTAB" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-HypreBiCGSTAB" );
    #ifdef AMP_USE_PETSC
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscFGMRES" );
            //        files.emplace_back(
            //        "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscBiCGSTAB" );
    #endif
#endif

#ifdef AMP_USE_PETSC
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-PetscCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-PetscFGMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-PetscBiCGSTAB" );
#endif

#ifdef AMP_USE_TRILINOS_ML
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-CG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-IPCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-FCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-GMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-FGMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-BiCGSTAB" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-TFQMR" );
    #ifdef AMP_USE_PETSC
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-PetscCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-PetscFGMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-PetscBiCGSTAB" );
    #endif
#endif

#ifdef AMP_USE_TRILINOS_MUELU
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-MueLu" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-MueLu-CG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-MueLu-IPCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-MueLu-FCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-MueLu-GMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-MueLu-FGMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-MueLu-BiCGSTAB" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-MueLu-TFQMR" );
    #ifdef AMP_USE_PETSC
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-MueLu-PetscCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-MueLu-PetscFGMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-MueLu-PetscBiCGSTAB" );
    #endif
#endif
    }

    std::vector<std::string> backends;
    backends.emplace_back( "serial" );
#if ( defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS ) )
    backends.emplace_back( "kokkos" );
#endif
#ifdef USE_DEVICE
    backends.emplace_back( "hip_cuda" );
#endif

    {
        PROFILE( "DRIVER::main(test loop)" );
        for ( auto &file : files ) {
            for ( auto &backend : backends )
                linearThermalTest( &ut, file, backend );
        }
    }

    ut.report();

    // build unique profile name to avoid collisions
    std::ostringstream ss;
    ss << "testLinSolveRobin_r" << std::setw( 3 ) << std::setfill( '0' )
       << AMP::AMPManager::getCommWorld().getSize();

    PROFILE_SAVE( ss.str() );

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
