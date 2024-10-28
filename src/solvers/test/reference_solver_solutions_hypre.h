#include "AMP/solvers/SolverStrategy.h"

#include <array>
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
// the reference number of iterations, final residual L2 norm maximum
// each version corresponds to a different architecture
#if ( defined USE_CUDA )
std::map<std::string, std::tuple<std::array<int, 3>, double>> conv_map{
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG",
      std::make_tuple( std::array<int, 3>( { 14, 14, 14 } ), REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG",
      std::make_tuple( std::array<int, 3>( { 7, 7, 7 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG",
      std::make_tuple( std::array<int, 3>( { 16, 16, 16 } ), REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG-CG",
      std::make_tuple( std::array<int, 3>( { 9, 9, 9 } ), REF_SOLVE_REL_LARGE ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRES",
      std::make_tuple( std::array<int, 3>( { 8, 8, 8 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FGMRES",
      std::make_tuple( std::array<int, 3>( { 8, 8, 8 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-BiCGSTAB",
      std::make_tuple( std::array<int, 3>( { 4, 4, 4 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-TFQMR",
      std::make_tuple( std::array<int, 3>( { 5, 5, 5 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscFGMRES",
      std::make_tuple( std::array<int, 3>( { 8, 8, 8 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreCG",
      std::make_tuple( std::array<int, 3>( { 8, 8, 8 } ), REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreCG",
      std::make_tuple( std::array<int, 3>( { 5, 5, 5 } ), REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-HypreCG",
      std::make_tuple( std::array<int, 3>( { 25, 25, 25 } ), REF_SOLVE_ABS ) }
};
#elif ( defined USE_HIP )
std::map<std::string, std::tuple<std::array<int, 3>, double>> conv_map{
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG",
      std::make_tuple( std::array<int, 3>( { 14, 14, 14 } ), REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG",
      std::make_tuple( std::array<int, 3>( { 7, 7, 7 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG",
      std::make_tuple( std::array<int, 3>( { 16, 16, 16 } ), REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG-CG",
      std::make_tuple( std::array<int, 3>( { 9, 9, 9 } ), REF_SOLVE_REL_LARGE ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRES",
      std::make_tuple( std::array<int, 3>( { 8, 8, 8 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FGMRES",
      std::make_tuple( std::array<int, 3>( { 8, 8, 8 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-BiCGSTAB",
      std::make_tuple( std::array<int, 3>( { 4, 3, 3 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-TFQMR",
      std::make_tuple( std::array<int, 3>( { 5, 5, 5 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscFGMRES",
      std::make_tuple( std::array<int, 3>( { 8, 8, 8 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreCG",
      std::make_tuple( std::array<int, 3>( { 8, 8, 8 } ), REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreCG",
      std::make_tuple( std::array<int, 3>( { 5, 5, 5 } ), REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-HypreCG",
      std::make_tuple( std::array<int, 3>( { 25, 25, 25 } ), REF_SOLVE_ABS ) }
};
#else
std::map<std::string, std::tuple<std::array<int, 3>, double>> conv_map{
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG",
      std::make_tuple( std::array<int, 3>( { 14, 14, 14 } ), REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG",
      std::make_tuple( std::array<int, 3>( { 7, 7, 7 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG",
      std::make_tuple( std::array<int, 3>( { 16, 16, 16 } ), REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG-CG",
      std::make_tuple( std::array<int, 3>( { 9, 9, 9 } ), REF_SOLVE_REL_LARGE ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRES",
      std::make_tuple( std::array<int, 3>( { 8, 8, 8 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FGMRES",
      std::make_tuple( std::array<int, 3>( { 8, 8, 8 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-BiCGSTAB",
      std::make_tuple( std::array<int, 3>( { 4, 4, 4 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-TFQMR",
      std::make_tuple( std::array<int, 3>( { 5, 5, 5 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscFGMRES",
      std::make_tuple( std::array<int, 3>( { 8, 8, 8 } ), REF_SOLVE_REL_SMALL ) },
    { "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreCG",
      std::make_tuple( std::array<int, 3>( { 8, 8, 8 } ), REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreCG",
      std::make_tuple( std::array<int, 3>( { 5, 5, 5 } ), REF_SOLVE_ABS ) },
    { "input_testLinearSolvers-LinearThermalRobin-HypreCG",
      std::make_tuple( std::array<int, 3>( { 25, 25, 25 } ), REF_SOLVE_ABS ) }
};
#endif

// Function to get the "solution" convergence rate and iteration count for the
// given input
static std::tuple<std::array<int, 3>, double> get_regression_solution( const std::string &input )
{
    auto it = conv_map.find( input );
    return ( it != conv_map.end() ) ?
               it->second :
               std::tuple<std::array<int, 3>, double>( std::array<int, 3>( { 0, 0, 0 } ), 0.0 );
}

static bool known_solution( const std::string &input )
{
    return ( conv_map.find( input ) != conv_map.end() );
}

static void checkConvergence( AMP::Solver::SolverStrategy *solver,
                              const std::string &inputFile,
                              const int comm_size,
                              AMP::UnitTest &ut )
{
    // allowed comm sizes are 1, 2, and 4.
    // Other sizes not covered by the test suite directly
    const bool allowed_size = comm_size == 1 || comm_size == 2 || comm_size == 4;
    const int size_idx      = comm_size == 4 ? 2 : ( comm_size == 2 ? 1 : 0 );

    const auto residualNorm = static_cast<double>( solver->getResidualNorm() );

    // Change this to check convergence status is on abstol or reltol.
    // Need to update solvers to set those flags first
    auto tolerance = static_cast<double>( solver->getAbsoluteTolerance() );

    if ( allowed_size && known_solution( inputFile ) ) {
        auto solution         = get_regression_solution( inputFile );
        const auto &ref_iters = std::get<0>( solution );
        const auto ref_iter   = ref_iters[size_idx];
        tolerance             = std::get<1>( solution );
        if ( ref_iter > 0 ) {
            const int iter = solver->getIterations();
            if ( iter != ref_iter || residualNorm > tolerance ) {
                AMP::pout << "FAILED: test_CellPreconditioner " << inputFile << std::endl;
                AMP::pout << "Iterations: " << iter << ", expected: " << ref_iter << std::endl;
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
