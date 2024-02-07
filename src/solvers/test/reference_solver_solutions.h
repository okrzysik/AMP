#include "AMP/solvers/SolverStrategy.h"

#include <iomanip>
#include <map>
#include <tuple>

// global data structure mapping an input file to a tuple of
// the reference number of iterations, final residual L2 norm,
// and the tolerance to which the final residual must match
std::map<std::string, std::tuple<int, double, double>> conv_map{
    { "input_testLinearSolvers-LinearThermalRobin-CG",
      std::make_tuple( 24, 1.68504e-10, 2.0e-10 ) },
    { "input_testLinearSolvers-LinearThermalRobin-GMRES",
      std::make_tuple( 25, 3.04412e-12, 4.0e-12 ) },
    { "input_testLinearSolvers-LinearThermalRobin-FGMRES",
      std::make_tuple( 25, 3.04412e-12, 4.0e-12 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BiCGSTAB",
      std::make_tuple( 18, 3.686753e-09, 1.0e-09 ) },
    { "input_testLinearSolvers-LinearThermalRobin-TFQMR",
      std::make_tuple( 20, 6.93865e-10, 1.0e-10 ) },
    { "input_testLinearSolvers-LinearThermalRobin-PetscFGMRES",
      std::make_tuple( 25, 6.62678786883557e-12, 6.0e-12 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG",
      std::make_tuple( 14, 3.90003330750106e-13, 4.0e-13 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG",
      std::make_tuple( 7, 2.04729351595174e-10, 3.0e-10 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRES",
      std::make_tuple( 8, 9.43856e-11, 1.0e-10 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FGMRES",
      std::make_tuple( 8, 9.47296e-11, 1.0e-10 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-BiCGSTAB",
      std::make_tuple( 4, 3.31261257647315e-12, 3.0e-12 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-TFQMR",
      std::make_tuple( 6, 7.93552e-10, 1.0e-11 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscFGMRES",
      std::make_tuple( 8, 9.414615740441e-11, 1.0e-10 ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreCG",
      std::make_tuple( 8, 1.58595265149432e-14, 2.0e-14 ) },
    { "input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreCG",
      std::make_tuple( 5, 1.20157627947399e-16, 1.0e-15 ) },
    { "input_testLinearSolvers-LinearThermalRobin-HypreCG",
      std::make_tuple( 25, 1.45084840704152e-15, 1.0e-14 ) },
    //    { "input_testLinearSolvers-LinearThermalRobin-ML",
    //      std::make_tuple( 14, 3.90003330750106e-13, 4.0e-13 ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-CG",
      std::make_tuple( 0, 2.55995674768853e-09, 2.0e-09 ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-GMRES",
      std::make_tuple( 1, 2.01677528883695e-09, 2.0e-9 ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-FGMRES",
      std::make_tuple( 1, 2.01677528883695e-09, 2.0e-9 ) },
    { "input_testLinearSolvers-LinearThermalRobin-ML-BiCGSTAB",
      std::make_tuple( 0, 1.03221e-21, 1.0e-12 ) },
    //    { "input_testLinearSolvers-LinearThermalRobin-ML-TFQMR",
    //      std::make_tuple( 6, 7.93552e-10, 1.0e-11 ) },
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
    const auto tolerance    = solver->getAbsoluteTolerance();

    if ( known_solution( inputFile ) ) {
        // Check the convergence rate to see if it changed
        auto solution = get_regression_solution( inputFile );
        auto ref_iter = std::get<0>( solution );
        auto ref_norm = std::get<1>( solution );
        auto ref_tol  = std::get<2>( solution );
        if ( ref_iter > 0 ) {
            int iter = solver->getIterations();
            if ( iter != ref_iter ||
                 fabs( static_cast<double>( residualNorm - ref_norm ) ) > ref_tol ) {
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
