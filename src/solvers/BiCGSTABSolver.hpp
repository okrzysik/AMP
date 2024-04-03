#include "AMP/operators/LinearOperator.h"
#include "AMP/solvers/BiCGSTABSolver.h"
#include "AMP/solvers/SolverFactory.h"

#include "ProfilerApp.h"

#include <cmath>
#include <limits>

namespace AMP::Solver {

/****************************************************************
 *  Constructors                                                 *
 ****************************************************************/
template<typename T>
BiCGSTABSolver<T>::BiCGSTABSolver() : d_restarts( 0 )
{
}

template<typename T>
BiCGSTABSolver<T>::BiCGSTABSolver( std::shared_ptr<SolverStrategyParameters> parameters )
    : SolverStrategy( parameters ), d_restarts( 0 )
{
    AMP_ASSERT( parameters );
    initialize( parameters );
}


/****************************************************************
 *  Destructor                                                   *
 ****************************************************************/
template<typename T>
BiCGSTABSolver<T>::~BiCGSTABSolver() = default;

/****************************************************************
 *  Initialize                                                   *
 ****************************************************************/
template<typename T>
void BiCGSTABSolver<T>::initialize( std::shared_ptr<const SolverStrategyParameters> parameters )
{
    AMP_ASSERT( parameters );
    auto db = parameters->d_db;
    getFromInput( db );

    if ( parameters->d_pNestedSolver ) {
        d_pPreconditioner = parameters->d_pNestedSolver;
    } else {
        if ( d_bUsesPreconditioner ) {
            auto pcName = db->getWithDefault<std::string>( "pc_solver_name", "Preconditioner" );
            std::shared_ptr<AMP::Database> outerDB;
            outerDB = db->keyExists( pcName ) ? db : parameters->d_global_db;
            AMP_INSIST( outerDB, "Outer database containing preconditioner is NULL" );
            auto pcDB       = outerDB->getDatabase( "Preconditioner" );
            auto parameters = std::make_shared<AMP::Solver::SolverStrategyParameters>( pcDB );
            parameters->d_pOperator = d_pOperator;
            d_pPreconditioner       = AMP::Solver::SolverFactory::create( parameters );
            AMP_ASSERT( d_pPreconditioner );
        }
    }
}

// Function to get values from input
template<typename T>
void BiCGSTABSolver<T>::getFromInput( std::shared_ptr<AMP::Database> db )
{
    d_bUsesPreconditioner = db->getWithDefault<bool>( "uses_preconditioner", false );
}

/****************************************************************
 *  Solve                                                        *
 * TODO: store convergence history, iterations, convergence reason
 ****************************************************************/
template<typename T>
void BiCGSTABSolver<T>::apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                               std::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE_START( "solve" );

    // Check input vector states
    AMP_ASSERT( ( f->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::UNCHANGED ) ||
                ( f->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::LOCAL_CHANGED ) );
    AMP_ASSERT( ( u->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::UNCHANGED ) ||
                ( u->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::LOCAL_CHANGED ) );

    // compute the norm of the rhs in order to compute
    // the termination criterion
    auto f_norm = static_cast<T>( f->L2Norm() );

    // if the rhs is zero we try to converge to the relative convergence
    if ( f_norm == static_cast<T>( 0.0 ) ) {
        f_norm = static_cast<T>( 1.0 );
    }

    const T terminate_tol = std::max( static_cast<T>( d_dRelativeTolerance * f_norm ),
                                      static_cast<T>( d_dAbsoluteTolerance ) );

    if ( d_iDebugPrintInfoLevel > 2 ) {
        std::cout << "BiCGSTABSolver<T>::solve: initial L2Norm of solution vector: " << u->L2Norm()
                  << std::endl;
        std::cout << "BiCGSTABSolver<T>::solve: initial L2Norm of rhs vector: " << f_norm
                  << std::endl;
    }

    if ( d_pOperator ) {
        registerOperator( d_pOperator );
    }

    // residual vector
    AMP::LinearAlgebra::Vector::shared_ptr res = f->clone();

    // compute the initial residual
    if ( d_bUseZeroInitialGuess ) {
        res->copyVector( f );
    } else {
        d_pOperator->residual( f, u, res );
    }

    // compute the current residual norm
    auto res_norm     = static_cast<T>( res->L2Norm() );
    auto r_tilde_norm = res_norm;

    if ( d_iDebugPrintInfoLevel > 0 ) {
        std::cout << "BiCGSTAB: initial residual " << res_norm << std::endl;
    }

    // return if the residual is already low enough
    if ( res_norm < terminate_tol ) {
        // provide a convergence reason
        // provide history (iterations, conv history etc)
        if ( d_iDebugPrintInfoLevel > 0 ) {
            std::cout << "BiCGSTABSolver<T>::solve: initial residual norm " << res_norm
                      << " is below convergence tolerance: " << terminate_tol << std::endl;
        }
        return;
    }

    // parameters in BiCGSTAB
    T alpha = static_cast<T>( 1.0 );
    T beta  = static_cast<T>( 0.0 );
    T omega = static_cast<T>( 1.0 );
    std::vector<T> rho( 2, static_cast<T>( 1.0 ) );
    NULL_USE( beta );

    // r_tilde is a non-zero initial direction chosen to be r
    std::shared_ptr<AMP::LinearAlgebra::Vector> r_tilde;
    // traditional choice is the initial residual
    r_tilde = res->clone();
    r_tilde->copyVector( res );

    auto p = res->clone();
    auto v = res->clone();
    p->zero();
    v->zero();

    std::shared_ptr<AMP::LinearAlgebra::Vector> p_hat, s, s_hat, t;

    for ( d_iNumberIterations = 0; d_iNumberIterations < d_iMaxIterations; ++d_iNumberIterations ) {

        rho[1] = static_cast<T>( r_tilde->dot( *res ) );

        auto angle = std::sqrt( std::fabs( rho[1] ) );
        auto eps   = std::numeric_limits<T>::epsilon();

        if ( angle < eps * r_tilde_norm ) {
            // the method breaks down as the vectors are orthogonal to r0
            // attempt to restart with a new r0
            u->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
            d_pOperator->residual( f, u, res );
            r_tilde->copyVector( res );
            p->copyVector( res );
            res_norm = static_cast<T>( res->L2Norm() );
            rho[1] = r_tilde_norm = res_norm;
            NULL_USE( rho[1] );
            d_restarts++;
            continue;
        }

        if ( d_iNumberIterations == 0 ) {

            p->copyVector( res );
        } else {

            beta = ( rho[1] / rho[0] ) * ( alpha / omega );
            p->axpy( -omega, *v, *p );
            p->axpy( beta, *p, *res );
        }

        if ( !p_hat ) {
            p_hat = u->clone();
            p_hat->zero();
        }

        // apply the preconditioner if it exists
        if ( d_bUsesPreconditioner ) {
            d_pPreconditioner->apply( p, p_hat );
        } else {
            p_hat->copyVector( p );
        }

        p_hat->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
        d_pOperator->apply( p_hat, v );

        alpha = static_cast<T>( r_tilde->dot( *v ) );
        AMP_ASSERT( alpha != static_cast<T>( 0.0 ) );
        alpha = rho[1] / alpha;

        if ( !s ) {
            s = res->clone();
        }
        s->axpy( -alpha, *v, *res );

        const auto s_norm = s->L2Norm();

        if ( s_norm < d_dRelativeTolerance ) {
            // early convergence
            u->axpy( alpha, *p_hat, *u );
            // add in code for final relative residual
            break;
        }

        if ( !s_hat ) {
            s_hat = u->clone();
            s_hat->zero();
        }

        // apply the preconditioner if it exists
        if ( d_bUsesPreconditioner ) {
            d_pPreconditioner->apply( s, s_hat );
        } else {
            s_hat->copyVector( s );
        }


        if ( !t ) {
            t = res->clone();
        }

        s_hat->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
        d_pOperator->apply( s_hat, t );

        auto t_sqnorm = static_cast<T>( t->dot( *t ) );
        auto t_dot_s  = static_cast<T>( t->dot( *s ) );
        omega = ( t_sqnorm == static_cast<T>( 0.0 ) ) ? static_cast<T>( 0.0 ) : t_dot_s / t_sqnorm;

        u->axpy( alpha, *p_hat, *u );
        u->axpy( omega, *s_hat, *u );

        res->axpy( -omega, *t, *s );

        // compute the current residual norm
        res_norm = static_cast<T>( res->L2Norm() );

        if ( d_iDebugPrintInfoLevel > 0 ) {
            std::cout << "BiCGSTAB: iteration " << ( d_iNumberIterations + 1 ) << ", residual "
                      << res_norm << std::endl;
        }

        // break if the residual is already low enough
        if ( res_norm < terminate_tol ) {
            // provide a convergence reason
            // provide history (iterations, conv history etc)
            break;
        }

        if ( omega == static_cast<T>( 0.0 ) ) {
            // this is a breakdown of the iteration
            // need to flag
            if ( d_iDebugPrintInfoLevel > 0 ) {
                std::cout << "BiCGSTABSolver<T>::solve: breakdown encountered, omega == 0"
                          << std::endl;
            }
            break;
        }

        rho[0] = rho[1];
    }

    u->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );

    if ( d_bComputeResidual ) {
        d_pOperator->residual( f, u, res );
        d_dResidualNorm = static_cast<T>( res->L2Norm() );
    } else
        d_dResidualNorm = res_norm;

    if ( d_iDebugPrintInfoLevel > 2 ) {
        std::cout << "L2Norm of solution: " << u->L2Norm() << std::endl;
    }

    PROFILE_STOP( "solve" );
}

/****************************************************************
 *  Function to set the register the operator                    *
 ****************************************************************/
template<typename T>
void BiCGSTABSolver<T>::registerOperator( std::shared_ptr<AMP::Operator::Operator> op )
{
    AMP_ASSERT( op );
    d_pOperator = op;
}

template<typename T>
void BiCGSTABSolver<T>::resetOperator(
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
} // namespace AMP::Solver
