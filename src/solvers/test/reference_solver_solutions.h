#include "AMP/solvers/SolverStrategy.h"

#include <iomanip>
#include <map>
#include <tuple>

// global data structure mapping an input file to a tuple containing
//  <0> reference number of iterations
//  <1> reference residual L2 norm, pass negative to compute internally
//  <2> testing tolerance on final L2
//  <3> strict checking
// The strict checking flag makes the pass/fail check behave in two ways:
//   Strict true:  Iterations to convergence must match reference exactly,
//                 and residual L2 norm must match reference within
//                 tolerance.
//   Strict false: Reference iterations and residual norm are treated as
//                 maximum allowed values. Testing tolerance unused.
// In either case the solver must report that it converged with reason
// ConvergedOnRelTol or ConvergedOnAbsTol.
std::map<std::string, std::tuple<int, double, double, bool>> conv_map{
    { "input_testLinearSolvers-LinearThermalRobin-CG",
      std::make_tuple( 24, 1.68504e-10, 3.0e-10, true ) },
    { "input_testLinearSolvers-LinearThermalRobin-IPCG",
      std::make_tuple( 24, 2.0e-10, 3.0e-10, true ) },
    { "input_testLinearSolvers-LinearThermalRobin-FCG",
      std::make_tuple( 24, 1.60876e-10, 5.0e-10, true ) },
    { "input_testLinearSolvers-LinearThermalRobin-CylMesh-CG",
      std::make_tuple( 87, 2.4806061595639e-09, 3.0e-10, true ) },
    { "input_testLinearSolvers-LinearThermalRobin-GMRES",
      std::make_tuple( 24, 3.04412e-12, 5.0e-12, true ) },
    { "input_testLinearSolvers-LinearThermalRobin-FGMRES",
      std::make_tuple( 24, 3.04412e-12, 5.0e-12, true ) },
    { "input_testLinearSolvers-LinearThermalRobin-BiCGSTAB",
      std::make_tuple( 18, 3.686753e-09, 1.0e-09, true ) },
    { "input_testLinearSolvers-LinearThermalRobin-TFQMR",
      std::make_tuple( 20, 6.0e-09, 5.0e-09, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-PetscFGMRES",
      std::make_tuple( 25, 6.62678786883557e-12, 6.0e-12, true ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML",
      std::make_tuple( 0, 1.1368884998304e-12, 2.0e-12, true ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-CG",
      std::make_tuple( 5, 1.27898867439474e-09, 2.0e-10, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-IPCG",
      std::make_tuple( 5, 1.279016778102248e-09, 2.0e-09, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-GMRES",
      std::make_tuple( 5, 1.220976553256886e-09, 2.0e-10, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-FGMRES",
      std::make_tuple( 5, 1.220976553256886e-09, 2.0e-10, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-BiCGSTAB",
      std::make_tuple( 3, 2.41617681713843e-09, 2.0e-09, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-TFQMR",
      std::make_tuple( 4, 1.89450365588275e-10, 2.0e-10, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-PetscFGMRES",
      std::make_tuple( 6, 1.217147923031935e-09, 2.0e-09, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-MueLu",
      std::make_tuple( 0, 3.947153362581674e-09, 2.0e-10, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-MueLu-CG",
      std::make_tuple( 27, 3.61364e-09, 2.0e-09, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-MueLu-IPCG",
      std::make_tuple( 9, 5.092835569229348e-10, 2.0e-10, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-MueLu-GMRES",
      std::make_tuple( 10, 2.099544001351115e-09, 2.0e-09, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-MueLu-FGMRES",
      std::make_tuple( 10, 2.099544001351115e-09, 2.0e-09, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-MueLu-BiCGSTAB",
      std::make_tuple( 4, 2.696506029673837e-09, 2.0e-09, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-MueLu-TFQMR",
      std::make_tuple( 6, 7.02626e-10, 2.0e-10, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-MueLu-PetscFGMRES",
      std::make_tuple( 10, 2.085761625746e-09, 2.0e-09, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG",
      std::make_tuple( 15, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG",
      std::make_tuple( 7, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-IPCG",
      std::make_tuple( 7, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FCG",
      std::make_tuple( 7, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG-FCG",
      std::make_tuple( 4, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG",
      std::make_tuple( 28, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG-CG",
      std::make_tuple( 11, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRES",
      std::make_tuple( 8, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-TFQMR",
      std::make_tuple( 5, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-BiCGSTAB",
      std::make_tuple( 4, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-GCR",
      std::make_tuple( 24, -1.0, 0.0, false ) },
    //    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-GMRES", (8, -1.0, 0.0,
    //    false) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FGMRES",
      std::make_tuple( 8, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-BiCGSTAB",
      std::make_tuple( 4, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-TFQMR",
      std::make_tuple( 6, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscFGMRES",
      std::make_tuple( 8, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreCG",
      std::make_tuple( 8, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreCG",
      std::make_tuple( 5, -1.0, 0.0, false ) },
    { "input_testLinearSolvers-LinearThermalRobin-HypreCG",
      std::make_tuple( 25, -1.0, 0.0, false ) }
};

// Function to get the "solution" convergence rate and iteration count for the
// given input
static std::tuple<int, double, double, bool> get_regression_solution( const std::string &input )
{
    auto it = conv_map.find( input );
    return ( it != conv_map.end() ) ? it->second :
                                      std::tuple<int, double, double, bool>( 0, 0.0, 0.0, false );
}

static bool known_solution( const std::string &input )
{
    return ( conv_map.find( input ) != conv_map.end() );
}

// Test for validity of solution by referencing map above and following
// the rules laid out there.
// Unrecognized input files just check if convergence reason is ok
static void checkConvergence( AMP::Solver::SolverStrategy *solver,
                              const std::string &inputFile,
                              AMP::UnitTest &ut )
{
    const auto final_norm = solver->getResidualNorm();
    const int iter        = solver->getIterations();
    const auto calc_norm  = std::fmax( static_cast<double>( solver->getAbsoluteTolerance() ),
                                      static_cast<double>( solver->getInitialResidual() ) *
                                          static_cast<double>( solver->getRelativeTolerance() ) );

    // Accept solution only if converged on absolute or relative tolerance
    const auto convReason = solver->getConvergenceStatus();
    const bool accept =
        convReason == AMP::Solver::SolverStrategy::SolverStatus::ConvergedOnRelTol ||
        convReason == AMP::Solver::SolverStrategy::SolverStatus::ConvergedOnAbsTol;

    if ( known_solution( inputFile ) ) {
        // Get the reference information
        auto [ref_iter, ref_norm, ref_tol, strict] = get_regression_solution( inputFile );
        // override reference norm if needed
        ref_norm = ref_norm > 0.0 ? ref_norm : calc_norm;
        // set pass/fail for each condition
        const bool pass_iters = strict ? iter == ref_iter : iter <= ref_iter;
        const auto norm_diff  = static_cast<double>( final_norm - ref_norm );
        const bool pass_norm  = strict ? std::fabs( norm_diff ) <= ref_tol : norm_diff <= 0.0;
        // Report passing or dump out information and fail
        if ( pass_iters && pass_norm && accept ) {
            ut.passes( "Passes convergence rate test" );
        } else if ( strict ) {
            AMP::pout << "FAILED: " << inputFile << std::endl;
            AMP::pout << "Iterations: " << iter << ", residual norm: " << std::setprecision( 15 )
                      << final_norm << std::endl;
            AMP::pout << "Expected: Iterations: " << ref_iter
                      << ", residual norm: " << std::setprecision( 15 ) << ref_norm << std::endl;
            AMP::pout << "Difference ( computed - reference ), iterations: " << ( iter - ref_iter )
                      << ", L2 norms: " << norm_diff << ", tolerance: " << ref_tol << std::endl;
            AMP::pout << "  Solver finished with status: " << solver->getConvergenceStatusString()
                      << std::endl;
            ut.failure( "FAILED: convergence rate test" );
        } else {
            AMP::pout << "FAILED: " << inputFile << std::endl;
            AMP::pout << "Iterations: " << iter << ", residual norm: " << std::setprecision( 15 )
                      << final_norm << std::endl;
            AMP::pout << "Maximum: Iterations: " << ref_iter
                      << ", residual norm: " << std::setprecision( 15 ) << ref_norm << std::endl;
            AMP::pout << "  Solver finished with status: " << solver->getConvergenceStatusString()
                      << std::endl;
            ut.failure( "FAILED: convergence rate test" );
        }
    } else {
        if ( accept ) {
            ut.passes( "Solver has converged for " + inputFile );
        } else {
            AMP::pout << "Solver has NOT converged for " << inputFile << std::endl;
            AMP::pout << "  Solver finished with status: " << solver->getConvergenceStatusString()
                      << std::endl;
            ut.failure( "Solver has NOT converged." );
        }
    }
}
