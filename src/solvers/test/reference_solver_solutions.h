#include "AMP/solvers/SolverStrategy.h"

#include <iomanip>
#include <map>
#include <tuple>

// global data structure mapping an input file to a tuple of
// the reference number of iterations, final residual L2 norm,
// and the tolerance to which the final residual must match
std::map<std::string, std::tuple<int, double, double>> conv_map{
    { "input_testLinearSolvers-LinearThermalRobin-CG",
      std::make_tuple( 24, 1.68504e-10, 3.0e-10 ) },
    { "input_testLinearSolvers-LinearThermalRobin-IPCG", std::make_tuple( 24, 2.0e-10, 2.0e-10 ) },
    { "input_testLinearSolvers-LinearThermalRobin-FCG",
      std::make_tuple( 24, 1.60876e-10, 5.0e-10 ) },
    { "input_testLinearSolvers-LinearThermalRobin-CylMesh-CG",
      std::make_tuple( 87, 2.4806061595639e-09, 3.0e-10 ) },
    { "input_testLinearSolvers-LinearThermalRobin-GMRES",
      std::make_tuple( 24, 3.04412e-12, 5.0e-12 ) },
    { "input_testLinearSolvers-LinearThermalRobin-FGMRES",
      std::make_tuple( 24, 3.04412e-12, 5.0e-12 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BiCGSTAB",
      std::make_tuple( 18, 3.686753e-09, 1.0e-09 ) },
    { "input_testLinearSolvers-LinearThermalRobin-TFQMR",
      std::make_tuple( 20, 1.51873e-09, 2.0e-09 ) },
    { "input_testLinearSolvers-LinearThermalRobin-PetscFGMRES",
      std::make_tuple( 25, 6.62678786883557e-12, 6.0e-12 ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML",
      std::make_tuple( 14, 3.90003330750106e-13, 4.0e-13 ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-CG",
      std::make_tuple( 0, 2.55995674768853e-09, 2.0e-09 ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-IPCG",
      std::make_tuple( 0, 2.55995674768853e-09, 2.0e-09 ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-GMRES",
      std::make_tuple( 4, 2.56859864182836e-09, 7.0e-08 ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-FGMRES",
      std::make_tuple( 1, 2.01677528883695e-09, 2.0e-9 ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-BiCGSTAB",
      std::make_tuple( 0, 1.03221e-21, 1.0e-12 ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-TFQMR",
      std::make_tuple( 3, 1.03439e-10, 2.0e-10 ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-PetscFGMRES",
      std::make_tuple( 1, 2.00917711909309e-09, 2.0e-09 ) }
};

// Function to get the "solution" convergence rate and iteration count for the
// given input
static std::tuple<int, double, double> get_regression_solution( const std::string &input )
{
    auto it = conv_map.find( input );
    return ( it != conv_map.end() ) ? it->second : std::tuple<int, double, double>( 0, 0.0, 0.0 );
}

static bool known_solution( const std::string &input )
{
    return ( conv_map.find( input ) != conv_map.end() );
}

static void checkConvergence( AMP::Solver::SolverStrategy *solver,
                              const std::string &inputFile,
                              AMP::UnitTest &ut )
{
    const auto residualNorm = solver->getResidualNorm();

    // Accept solution only if converged on absolute or relative tolerance
    const auto convReason = solver->getConvergenceStatus();
    const bool accept =
        convReason == AMP::Solver::SolverStrategy::SolverStatus::ConvergedOnRelTol ||
        convReason == AMP::Solver::SolverStrategy::SolverStatus::ConvergedOnAbsTol;

    if ( known_solution( inputFile ) ) {
        // Check the convergence rate to see if it changed
        auto solution = get_regression_solution( inputFile );
        auto ref_iter = std::get<0>( solution );
        auto ref_norm = std::get<1>( solution );
        auto ref_tol  = std::get<2>( solution );
        if ( ref_iter > 0 ) {
            int iter = solver->getIterations();
            if ( iter != ref_iter ||
                 fabs( static_cast<double>( residualNorm - ref_norm ) ) > ref_tol || !accept ) {
                AMP::pout << "FAILED: test_CellPreconditioner " << inputFile << std::endl;
                AMP::pout << "Iterations: " << iter
                          << ", residual norm: " << std::setprecision( 15 ) << residualNorm
                          << std::endl;
                AMP::pout << "Expected: Iterations: " << ref_iter
                          << ", residual norm: " << std::setprecision( 15 ) << ref_norm
                          << std::endl;
                AMP::pout << "Difference ( computed - reference ), iterations: "
                          << ( iter - ref_iter ) << ", L2 norms: " << ( residualNorm - ref_norm )
                          << ", tolerance: " << ref_tol << std::endl;
                AMP::pout << "  Solver finished with status: "
                          << solver->getConvergenceStatusString() << std::endl;
                ut.failure( "FAILED: convergence rate test" );
            } else {
                ut.passes( "Passes convergence rate test" );
            }
        }
    } else {
        if ( accept ) {
            ut.passes( "Solver has converged." );
        } else {
            AMP::pout << "Solver has NOT converged." << std::endl;
            AMP::pout << "  Solver finished with status: " << solver->getConvergenceStatusString()
                      << std::endl;
            ut.failure( "Solver has NOT converged." );
        }
    }
}
