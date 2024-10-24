#include "AMP/solvers/SolverStrategy.h"

#include <iomanip>
#include <map>
#include <tuple>

// All test cases solve one of two problems
// The relative tolerance is scaled by initial residual norm
// but solvers don't always report their convergence status
//   => three possible criteria for convergence:
//         - Absolute tolerance, always 1.0e-12 rarely triggered
//         - Relative tolerance small prob, always 1.0e-12 * initial residual
//         - Relative tolerance large prob, always 1.0e-12 * initial residual
// Note: small refers to system size, not the residual size
const double REF_SOLVE_ABS       = 1.0e-12;
const double REF_SOLVE_REL_SMALL = 1.0e-12 * 4419.95;
const double REF_SOLVE_REL_LARGE = 1.0e-12 * 3092.94;

// global data structure mapping an input file to a tuple of
// the reference number of iterations, final residual L2 norm,
// and the tolerance to which the final residual must match
std::map<std::string, std::tuple<int, double>> conv_map{
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG",
      std::make_tuple( 14, REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG",
      std::make_tuple( 7, REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG",
      std::make_tuple( 16, REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG-CG",
      std::make_tuple( 9, REF_SOLVE_REL_LARGE ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRES",
      std::make_tuple( 8, REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FGMRES",
      std::make_tuple( 8, REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-BiCGSTAB",
      std::make_tuple( 4, REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-TFQMR",
      std::make_tuple( 5, REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscFGMRES",
      std::make_tuple( 8, REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreCG",
      std::make_tuple( 8, REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreCG",
      std::make_tuple( 5, REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-HypreCG", std::make_tuple( 25, REF_SOLVE_ABS ) }
};

// Function to get the "solution" convergence rate and iteration count for the
// given input
static std::tuple<int, double> get_regression_solution( const std::string &input )
{
    auto it = conv_map.find( input );
    return ( it != conv_map.end() ) ? it->second : std::tuple<int, double>( 0, 0.0 );
}

static bool known_solution( const std::string &input )
{
    return ( conv_map.find( input ) != conv_map.end() );
}

static void checkConvergence( AMP::Solver::SolverStrategy *solver,
                              const std::string &inputFile,
                              AMP::UnitTest &ut )
{
    const auto residualNorm = static_cast<double>( solver->getResidualNorm() );

    // Change this to check convergence status is on abstol or reltol.
    // Need to update solvers to set those flags first
    auto tolerance = static_cast<double>( solver->getAbsoluteTolerance() );

    if ( known_solution( inputFile ) ) {
        auto solution       = get_regression_solution( inputFile );
        const auto ref_iter = std::get<0>( solution );
        tolerance           = std::get<1>( solution );
        if ( ref_iter > 0 ) {
            const int iter = solver->getIterations();
            if ( iter > ref_iter || residualNorm > tolerance ) {
                AMP::pout << "FAILED: test_CellPreconditioner " << inputFile << std::endl;
                AMP::pout << "Iterations: " << iter << ", max allowed: " << ref_iter << std::endl;
                AMP::pout << "Residual norm: " << residualNorm << ", max allowed: " << tolerance
                          << std::endl;
                ut.failure( "FAILED: convergence rate test" );
            } else {
                ut.passes( "Passes convergence rate test" );
            }
        }
    } else {
        if ( residualNorm < tolerance ) {
            ut.passes( "Solver has converged." );
        } else {
            AMP::pout << "Solver has NOT converged." << std::endl;
            AMP::pout << "Residual norm: " << residualNorm << ", max allowed: " << tolerance
                      << std::endl;
            ut.failure( "Solver has NOT converged." );
        }
    }
}
