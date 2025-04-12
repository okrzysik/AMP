#include "AMP/solvers/SolverStrategy.h"

#include <iomanip>
#include <map>

// global data structure mapping an input file to a tuple of
// the reference number of iterations, final residual L2 norm maximum
// each version corresponds to a different architecture
std::map<std::string, int> conv_map{
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG", 15 },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG", 7 },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-IPCG", 7 },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FCG", 7 },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG-FCG", 4 },
    { "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG",
      28 }, // max is 18 looking just at GPU runs...
    { "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG-CG", 11 },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRES", 8 },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-TFQMR", 5 },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-BiCGSTAB", 4 },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-GCR", 24 },
    //    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-GMRES", 8 },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FGMRES", 8 },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-BiCGSTAB", 4 },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-TFQMR", 6 },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscFGMRES", 8 },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreCG", 8 },
    { "input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreCG", 5 },
    { "input_testLinearSolvers-LinearThermalRobin-HypreCG", 25 }
};

// Function to get the "solution" convergence rate and iteration count for the
// given input
static int getReferenceIterations( const std::string &input )
{
    auto it = conv_map.find( input );
    return ( it != conv_map.end() ) ? it->second : 0;
}

static void checkConvergence( AMP::Solver::SolverStrategy *solver,
                              const std::string &inputFile,
                              AMP::UnitTest &ut )
{
    // Get bound on allowable residual norm, just for reporting info on failure
    const auto initRes = static_cast<double>( solver->getInitialResidual() );
    const auto absTol  = static_cast<double>( solver->getAbsoluteTolerance() );
    const auto relTol  = static_cast<double>( solver->getRelativeTolerance() );
    const auto stopTol = std::fmax( absTol, initRes * relTol );

    // get final residual norm
    const auto residualNorm = static_cast<double>( solver->getResidualNorm() );

    // Accept solution only if converged on absolute or relative tolerance
    const auto convReason = solver->getConvergenceStatus();
    const bool accept =
        convReason == AMP::Solver::SolverStrategy::SolverStatus::ConvergedOnRelTol ||
        convReason == AMP::Solver::SolverStrategy::SolverStatus::ConvergedOnAbsTol;

    const auto refIter = getReferenceIterations( inputFile );
    if ( refIter > 0 ) {
        const int iter = solver->getIterations();
        if ( iter > refIter || !accept ) {
            AMP::pout << "FAILED: test_CellPreconditioner " << inputFile << std::endl;
            AMP::pout << "  Iterations: " << iter << ", max allowed: " << refIter << std::endl;
            AMP::pout << "  Residual norm: " << residualNorm << ", max allowed: " << stopTol
                      << std::endl;
            AMP::pout << "  Solver finished with status: " << solver->getConvergenceStatusString()
                      << std::endl;
            ut.failure( "FAILED: convergence rate test" );
        } else {
            ut.passes( "Passes convergence rate test" );
        }
    } else if ( accept ) {
        ut.passes( "Solver has converged." );
    } else {
        AMP::pout << "Solver has NOT converged." << std::endl;
        AMP::pout << "Residual norm: " << residualNorm << ", max allowed: " << stopTol << std::endl;
        AMP::pout << "  Solver finished with status: " << solver->getConvergenceStatusString()
                  << std::endl;
        ut.failure( "Solver has NOT converged." );
    }
}
