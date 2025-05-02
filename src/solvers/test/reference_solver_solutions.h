#include "AMP/solvers/SolverStrategy.h"

#include <iomanip>
#include <map>
#include <tuple>

// Function to get the "solution" convergence rate and iteration count for the
// given input
static std::tuple<int, double, double, bool>
get_regression_solution( std::shared_ptr<const AMP::Database> input_db )
{
    auto db        = input_db->getDatabase( "Reference" );
    auto its       = db->getWithDefault<int>( "iterations", 0 );
    auto res_norm  = db->getWithDefault<double>( "res_l2_norm", -1.0 );
    auto tolerance = db->getWithDefault<double>( "tolerance", 0.0 );
    auto strict    = db->getWithDefault<bool>( "strict", false );
    return std::tuple<int, double, double, bool>( its, res_norm, tolerance, strict );
}

// Test for validity of solution by referencing map above and following
// the rules laid out there.
// Unrecognized input files just check if convergence reason is ok
static void checkConvergence( AMP::Solver::SolverStrategy *solver,
                              std::shared_ptr<const AMP::Database> input_db,
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

    if ( input_db->keyExists( "Reference" ) ) {
        // Get the reference information
        auto [ref_iter, ref_norm, ref_tol, strict] = get_regression_solution( input_db );
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
