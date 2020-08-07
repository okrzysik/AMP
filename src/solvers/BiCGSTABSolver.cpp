#include "AMP/solvers/BiCGSTABSolver.h"
#include "AMP/operators/LinearOperator.h"
#include "AMP/solvers/KrylovSolverParameters.h"

#include "ProfilerApp.h"

#include <cmath>
#include <limits>

namespace AMP {
namespace Solver {

/****************************************************************
 *  Constructors                                                 *
 ****************************************************************/
BiCGSTABSolver::BiCGSTABSolver() : d_restarts( 0 ) {}

BiCGSTABSolver::BiCGSTABSolver( std::shared_ptr<SolverStrategyParameters> parameters )
    : SolverStrategy( parameters ), d_restarts( 0 )
{
    AMP_ASSERT( parameters.get() != nullptr );

    // Initialize
    initialize( parameters );
}


/****************************************************************
 *  Destructor                                                   *
 ****************************************************************/
BiCGSTABSolver::~BiCGSTABSolver() = default;

/****************************************************************
 *  Initialize                                                   *
 ****************************************************************/
void BiCGSTABSolver::initialize( std::shared_ptr<SolverStrategyParameters> const params )
{
    auto parameters = std::dynamic_pointer_cast<KrylovSolverParameters>( params );
    AMP_ASSERT( parameters.get() != nullptr );
    d_comm = parameters->d_comm;
    AMP_ASSERT( !d_comm.isNull() );

    d_pPreconditioner = parameters->d_pPreconditioner;

    getFromInput( parameters->d_db );

    if ( d_pOperator.get() != nullptr ) {
        registerOperator( d_pOperator );
    }
}

// Function to get values from input
void BiCGSTABSolver::getFromInput( std::shared_ptr<AMP::Database> db )
{

    d_dRelativeTolerance  = db->getWithDefault<double>( "relative_tolerance", 1.0e-9 );
    d_iMaxIterations      = db->getWithDefault<double>( "max_iterations", 1000 );
    d_bUsesPreconditioner = db->getWithDefault( "uses_preconditioner", false );
}

/****************************************************************
 *  Solve                                                        *
 * TODO: store convergence history, iterations, convergence reason
 ****************************************************************/
void BiCGSTABSolver::solve( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                            std::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE_START( "solve" );

    // Check input vector states
    AMP_ASSERT(
        ( f->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( f->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );
    AMP_ASSERT(
        ( u->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( u->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );

    // compute the norm of the rhs in order to compute
    // the termination criterion
    double f_norm = f->L2Norm( f );

    // if the rhs is zero we try to converge to the relative convergence
    if ( f_norm == 0.0 ) {
        f_norm = 1.0;
    }

    const double terminate_tol = d_dRelativeTolerance * f_norm;

    if ( d_iDebugPrintInfoLevel > 2 ) {
        std::cout << "BiCGSTABSolver::solve: initial L2Norm of solution vector: " << u->L2Norm( u )
                  << std::endl;
        std::cout << "BiCGSTABSolver::solve: initial L2Norm of rhs vector: " << f_norm << std::endl;
    }

    if ( d_pOperator.get() != nullptr ) {
        registerOperator( d_pOperator );
    }

    // residual vector
    AMP::LinearAlgebra::Vector::shared_ptr res = f->cloneVector();

    // compute the initial residual
    if ( d_bUseZeroInitialGuess ) {
        res->copyVector( f );
    } else {
        d_pOperator->residual( f, u, res );
    }

    // compute the current residual norm
    double res_norm     = res->L2Norm( res );
    double r_tilde_norm = res_norm;

    if ( d_iDebugPrintInfoLevel > 0 ) {
        std::cout << "BiCGSTAB: initial residual " << res_norm << std::endl;
    }

    // return if the residual is already low enough
    if ( res_norm < terminate_tol ) {
        // provide a convergence reason
        // provide history (iterations, conv history etc)
        if ( d_iDebugPrintInfoLevel > 0 ) {
            std::cout << "BiCGSTABSolver::solve: initial residual norm " << res_norm
                      << " is below convergence tolerance: " << terminate_tol << std::endl;
        }
        return;
    }

    // parameters in BiCGSTAB
    double alpha = 1.0;
    double beta  = 0.0;
    double omega = 1.0;
    std::vector<double> rho( 2, 1.0 );
    NULL_USE( beta );

    // r_tilde is a non-zero initial direction chosen to be r
    std::shared_ptr<AMP::LinearAlgebra::Vector> r_tilde;
    // traditional choice is the initial residual
    r_tilde = res->cloneVector();
    r_tilde->copyVector( res );

    auto p = res->cloneVector();
    auto v = res->cloneVector();
    p->zero( p );
    v->zero( v );

    for ( auto iter = 0; iter < d_iMaxIterations; ++iter ) {

        rho[1] = r_tilde->dot( res, r_tilde );

        double angle = std::sqrt( std::fabs( rho[1] ) );
        double eps   = std::numeric_limits<double>::epsilon();

        if ( angle < eps * r_tilde_norm ) {
            // the method breaks down as the vectors are orthogonal to r0
            // attempt to restart with a new r0
            d_pOperator->residual( f, u, res );
            r_tilde->copyVector( res );
            res_norm = res->L2Norm( res );
            rho[1] = r_tilde_norm = res_norm;
            NULL_USE( rho[1] );
            d_restarts++;
            break;
        }

        if ( iter == 1 ) {

            p->copyVector( res );
        } else {

            beta = ( rho[1] / rho[0] ) * ( alpha / omega );
            p->axpy( -omega, v, p, p );
            p->axpy( beta, p, res, p );
        }

        std::shared_ptr<AMP::LinearAlgebra::Vector> p_hat = u->cloneVector();

        // apply the preconditioner if it exists
        if ( d_bUsesPreconditioner ) {
            d_pPreconditioner->solve( p, p_hat );
        } else {
            p_hat->copyVector( p );
        }

        d_pOperator->apply( p_hat, v );

        alpha = r_tilde->dot( v, r_tilde );
        AMP_ASSERT( alpha != 0.0 );
        alpha = rho[1] / alpha;

        std::shared_ptr<AMP::LinearAlgebra::Vector> s = res->cloneVector();
        s->axpy( -alpha, v, res, s );

        const double s_norm = s->L2Norm( s );

        if ( s_norm < d_dRelativeTolerance ) {
            // early convergence
            u->axpy( alpha, p_hat, u, u );
            // add in code for final relative residual
            break;
        }

        std::shared_ptr<AMP::LinearAlgebra::Vector> s_hat = u->cloneVector();

        // apply the preconditioner if it exists
        if ( d_bUsesPreconditioner ) {
            d_pPreconditioner->solve( s, s_hat );
        } else {
            s_hat->copyVector( s );
        }

        std::shared_ptr<AMP::LinearAlgebra::Vector> t = res->cloneVector();
        d_pOperator->apply( s_hat, t );

        const double t_sqnorm = t->dot( t, t );
        omega                 = ( t_sqnorm == 0.0 ) ? 0.0 : t->dot( s, t ) / t_sqnorm;

        u->axpy( alpha, p_hat, u, u );
        u->axpy( omega, s_hat, u, u );

        res->axpy( -omega, t, s, res );

        // compute the current residual norm
        res_norm = res->L2Norm( res );

        if ( d_iDebugPrintInfoLevel > 0 ) {
            std::cout << "BiCGSTAB: iteration " << ( iter + 1 ) << ", residual " << res_norm
                      << std::endl;
        }

        // break if the residual is already low enough
        if ( res_norm < terminate_tol ) {
            // provide a convergence reason
            // provide history (iterations, conv history etc)
            break;
        }

        if ( omega == 0.0 ) {
            // this is a breakdown of the iteration
            // need to flag
            if ( d_iDebugPrintInfoLevel > 0 ) {
                std::cout << "BiCGSTABSolver::solve: breakdown encountered, omega == 0"
                          << std::endl;
            }
            break;
        }

        rho[0] = rho[1];
    }

    if ( d_iDebugPrintInfoLevel > 2 ) {
        std::cout << "L2Norm of solution: " << u->L2Norm( u ) << std::endl;
    }

    PROFILE_STOP( "solve" );
}

/****************************************************************
 *  Function to set the register the operator                    *
 ****************************************************************/
void BiCGSTABSolver::registerOperator( const std::shared_ptr<AMP::Operator::Operator> op )
{
    AMP_ASSERT( op.get() != nullptr );

    d_pOperator = op;

    std::shared_ptr<AMP::Operator::LinearOperator> linearOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearOperator>( op );
    AMP_ASSERT( linearOperator.get() != nullptr );
}
void BiCGSTABSolver::resetOperator(
    const std::shared_ptr<AMP::Operator::OperatorParameters> params )
{
    if ( d_pOperator.get() != nullptr ) {
        d_pOperator->reset( params );
    }

    // should add a mechanism for the linear operator to provide updated parameters for the
    // preconditioner operator
    // though it's unclear where this might be necessary
    if ( d_pPreconditioner.get() != nullptr ) {
        d_pPreconditioner->resetOperator( params );
    }
}
} // namespace Solver
} // namespace AMP
