#include "AMP/solvers/SolverStrategy.h"

#include <iomanip>
#include <map>
#include <tuple>

// global data structure mapping an input file to a tuple of
// the reference number of iterations, final residual L2 norm,
// and the tolerance to which the final residual must match
std::map<std::string, std::tuple<int, double, double>> conv_map{
    { "input_testLinearSolvers-LinearThermalRobin-CG",
      std::make_tuple( 25, 1.68504e-10, 1.0e-10 ) },
    { "input_testLinearSolvers-LinearThermalRobin-GMRES",
      std::make_tuple( 25, 5.68264e-12, 1.0e-12 ) },
    { "input_testLinearSolvers-LinearThermalRobin-FGMRES",
      std::make_tuple( 25, 5.68264e-12, 1.0e-12 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BiCGSTAB",
      std::make_tuple( 19, 3.68678e-09, 1.0e-09 ) },
    { "input_testLinearSolvers-LinearThermalRobin-TFQMR",
      std::make_tuple( 21, 6.93865e-10, 1.0e-10 ) },
    { "input_testLinearSolvers-LinearThermalRobin-PetscFGMRES",
      std::make_tuple( 25, 9.059409431777e-12, 1.0e-12 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG",
      std::make_tuple( 10, 4.332373580517e-06, 1.0e-06 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG",
      std::make_tuple( 8, 9.9041e-11, 1.0e-12 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRES",
      std::make_tuple( 8, 9.43856e-11, 1.0e-12 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FGMRES",
      std::make_tuple( 8, 9.47296e-11, 1.0e-12 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-BiCGSTAB",
      std::make_tuple( 5, 3.64001e-12, 1.0e-13 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-TFQMR",
      std::make_tuple( 6, 7.93552e-10, 1.0e-11 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreCG",
      std::make_tuple( 8, 9.88953e-11, 1.0e-12 ) },
    { "input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreCG",
      std::make_tuple( 5, 2.22094e-12, 1.0e-13 ) },
    { "input_testLinearSolvers-LinearThermalRobin-HypreCG",
      std::make_tuple( 20, 2.678101722924e-04, 1.0e-04 ) }
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
    const auto tolerance    = solver->getAbsoluteTolerance();

    if ( known_solution( inputFile ) ) {
        // Check the convergence rate to see if it changed
        auto solution = get_regression_solution( inputFile );
        auto ref_iter = std::get<0>( solution );
        auto ref_norm = std::get<1>( solution );
        auto ref_tol  = std::get<2>( solution );
        if ( ref_iter > 0 ) {
            int iter = solver->getIterations();
            if ( iter != ref_iter || fabs( residualNorm - ref_norm ) > ref_tol ) {
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
            AMP::pout << "Residual norm: " << residualNorm << " is not smaller than tolerance "
                      << tolerance << "." << std::endl;
            ut.failure( "Solver has NOT converged." );
        }
    }
}
