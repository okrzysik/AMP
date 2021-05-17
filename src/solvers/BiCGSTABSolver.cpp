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
    AMP_ASSERT( parameters );

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
void BiCGSTABSolver::initialize( std::shared_ptr<const SolverStrategyParameters> params )
{
    auto parameters = std::dynamic_pointer_cast<const KrylovSolverParameters>( params );
    AMP_ASSERT( parameters );
    d_comm = parameters->d_comm;
    AMP_ASSERT( !d_comm.isNull() );

    d_pPreconditioner = parameters->d_pPreconditioner;

    getFromInput( parameters->d_db );

    if ( d_pOperator ) {
        registerOperator( d_pOperator );
    }
}

// Function to get values from input
void BiCGSTABSolver::getFromInput( std::shared_ptr<AMP::Database> db )
{

    d_iMaxIterations      = db->getWithDefault<double>( "max_iterations", 1000 );
    d_bUsesPreconditioner = db->getWithDefault( "uses_preconditioner", false );
}

/****************************************************************
 *  Solve                                                        *
 * TODO: store convergence history, iterations, convergence reason
 ****************************************************************/
void BiCGSTABSolver::apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                            std::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE_START( "solve" );

    // Check input vector states
    AMP_ASSERT(
        ( f->getUpdateStatus() == AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED ) ||
        ( f->getUpdateStatus() == AMP::LinearAlgebra::VectorData::UpdateState::LOCAL_CHANGED ) );
    AMP_ASSERT(
        ( u->getUpdateStatus() == AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED ) ||
        ( u->getUpdateStatus() == AMP::LinearAlgebra::VectorData::UpdateState::LOCAL_CHANGED ) );

    // compute the norm of the rhs in order to compute
    // the termination criterion
    double f_norm = static_cast<double>( f->L2Norm() );

    // if the rhs is zero we try to converge to the relative convergence
    if ( f_norm == 0.0 ) {
        f_norm = 1.0;
    }

    const double terminate_tol = d_dRelativeTolerance * f_norm;

    if ( d_iDebugPrintInfoLevel > 2 ) {
        std::cout << "BiCGSTABSolver::solve: initial L2Norm of solution vector: " << u->L2Norm()
                  << std::endl;
        std::cout << "BiCGSTABSolver::solve: initial L2Norm of rhs vector: " << f_norm << std::endl;
    }

    if ( d_pOperator ) {
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
    double res_norm     = static_cast<double>( res->L2Norm() );
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
    p->zero();
    v->zero();

    for ( auto iter = 0; iter < d_iMaxIterations; ++iter ) {

        rho[1] = static_cast<double>( r_tilde->dot( *res ) );

        double angle = std::sqrt( std::fabs( rho[1] ) );
        double eps   = std::numeric_limits<double>::epsilon();

        if ( angle < eps * r_tilde_norm ) {
            // the method breaks down as the vectors are orthogonal to r0
            // attempt to restart with a new r0
            d_pOperator->residual( f, u, res );
            r_tilde->copyVector( res );
            res_norm = static_cast<double>( res->L2Norm() );
            rho[1] = r_tilde_norm = res_norm;
            NULL_USE( rho[1] );
            d_restarts++;
            break;
        }

        if ( iter == 1 ) {

            p->copyVector( res );
        } else {

            beta = ( rho[1] / rho[0] ) * ( alpha / omega );
            p->axpy( -omega, *v, *p );
            p->axpy( beta, *p, *res );
        }

        auto p_hat = u->cloneVector();

        // apply the preconditioner if it exists
        if ( d_bUsesPreconditioner ) {
            d_pPreconditioner->apply( p, p_hat );
        } else {
            p_hat->copyVector( p );
        }

        d_pOperator->apply( p_hat, v );

        alpha = static_cast<double>( r_tilde->dot( *v ) );
        AMP_ASSERT( alpha != 0.0 );
        alpha = rho[1] / alpha;

        auto s = res->cloneVector();
        s->axpy( -alpha, *v, *res );

        const double s_norm = static_cast<double>( s->L2Norm() );

        if ( s_norm < d_dRelativeTolerance ) {
            // early convergence
            u->axpy( alpha, *p_hat, *u );
            // add in code for final relative residual
            break;
        }

        auto s_hat = u->cloneVector();

        // apply the preconditioner if it exists
        if ( d_bUsesPreconditioner ) {
            d_pPreconditioner->apply( s, s_hat );
        } else {
            s_hat->copyVector( s );
        }

        auto t = res->cloneVector();
        d_pOperator->apply( s_hat, t );

        double t_sqnorm = static_cast<double>( t->dot( *t ) );
        double t_dot_s  = static_cast<double>( t->dot( *s ) );
        omega           = ( t_sqnorm == 0.0 ) ? 0.0 : t_dot_s / t_sqnorm;

        u->axpy( alpha, *p_hat, *u );
        u->axpy( omega, *s_hat, *u );

        res->axpy( -omega, *t, *s );

        // compute the current residual norm
        res_norm = static_cast<double>( res->L2Norm() );

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
        std::cout << "L2Norm of solution: " << u->L2Norm() << std::endl;
    }

    PROFILE_STOP( "solve" );
}

/****************************************************************
 *  Function to set the register the operator                    *
 ****************************************************************/
void BiCGSTABSolver::registerOperator( std::shared_ptr<AMP::Operator::Operator> op )
{
    AMP_ASSERT( op );

    d_pOperator = op;

    std::shared_ptr<AMP::Operator::LinearOperator> linearOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearOperator>( op );
    AMP_ASSERT( linearOperator );
}
void BiCGSTABSolver::resetOperator(
    std::shared_ptr<const AMP::Operator::OperatorParameters> params )
{
    if ( d_pOperator ) {
        d_pOperator->reset( params );
    }

    // should add a mechanism for the linear operator to provide updated parameters for the
    // preconditioner operator
    // though it's unclear where this might be necessary
    if ( d_pPreconditioner ) {
        d_pPreconditioner->resetOperator( params );
    }
}
} // namespace Solver
} // namespace AMP
