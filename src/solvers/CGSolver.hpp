#include "AMP/operators/LinearOperator.h"
#include "AMP/solvers/CGSolver.h"
#include "AMP/solvers/SolverFactory.h"
#include "ProfilerApp.h"

namespace AMP::Solver {

/****************************************************************
 *  Constructors                                                 *
 ****************************************************************/

template<typename T>
CGSolver<T>::CGSolver( std::shared_ptr<AMP::Solver::SolverStrategyParameters> parameters )
    : SolverStrategy( parameters )
{
    AMP_ASSERT( parameters );

    // Initialize
    initialize( parameters );
}

/****************************************************************
 *  Initialize                                                   *
 ****************************************************************/
template<typename T>
void CGSolver<T>::initialize(
    std::shared_ptr<const AMP::Solver::SolverStrategyParameters> parameters )
{
    AMP_ASSERT( parameters );

    auto db = parameters->d_db;
    getFromInput( db );

    if ( parameters->d_pNestedSolver ) {
        d_pPreconditioner = parameters->d_pNestedSolver;
    } else {
        if ( d_bUsesPreconditioner ) {
            auto pcName  = db->getWithDefault<std::string>( "pc_solver_name", "Preconditioner" );
            auto outerDB = db->keyExists( pcName ) ? db : parameters->d_global_db;
            if ( outerDB ) {
                auto pcDB       = outerDB->getDatabase( pcName );
                auto parameters = std::make_shared<AMP::Solver::SolverStrategyParameters>( pcDB );
                parameters->d_pOperator = d_pOperator;
                d_pPreconditioner       = AMP::Solver::SolverFactory::create( parameters );
                AMP_ASSERT( d_pPreconditioner );
            }
        }
    }
}

// Function to get values from input
template<typename T>
void CGSolver<T>::getFromInput( std::shared_ptr<AMP::Database> db )
{
    d_dDivergenceTolerance = db->getWithDefault<T>( "divergence_tolerance", 1.0e+03 );
    d_bUsesPreconditioner  = db->getWithDefault<bool>( "uses_preconditioner", false );
}

/****************************************************************
 *  Solve                                                        *
 * TODO: store convergence history, iterations, convergence reason
 ****************************************************************/
template<typename T>
void CGSolver<T>::apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                         std::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE( "solve" );

    // Check input vector states
    AMP_ASSERT( ( f->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::UNCHANGED ) ||
                ( f->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::LOCAL_CHANGED ) );
    AMP_ASSERT( ( u->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::UNCHANGED ) ||
                ( u->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::LOCAL_CHANGED ) );

    // compute the norm of the rhs in order to compute
    // the termination criterion
    const auto f_norm = static_cast<T>( f->L2Norm() );

    // enhance with convergence reason, number of iterations etc
    if ( f_norm == static_cast<T>( 0.0 ) ) {
        u->zero();
        d_ConvergenceStatus = SolverStatus::ConvergedOnAbsTol;
        return;
    }
    const auto terminate_tol = std::max( static_cast<T>( d_dRelativeTolerance * f_norm ),
                                         static_cast<T>( d_dAbsoluteTolerance ) );

    // z will store r when a preconditioner is not present
    // and will store the result of a preconditioner solve
    // when a preconditioner is present
    std::shared_ptr<AMP::LinearAlgebra::Vector> z;

    // residual vector
    AMP::LinearAlgebra::Vector::shared_ptr r = f->clone();

    // compute the initial residual
    if ( d_bUseZeroInitialGuess ) {
        u->zero();
        r->copyVector( f );
    } else {
        d_pOperator->residual( f, u, r );
    }

    if ( d_iDebugPrintInfoLevel > 1 ) {
        AMP::pout << "CGSolver<T>::solve: initial L2Norm of solution vector: " << u->L2Norm()
                  << std::endl;
        AMP::pout << "CGSolver<T>::solve: initial L2Norm of rhs vector: " << f_norm << std::endl;
    }
    // compute the current residual norm
    auto current_res = static_cast<T>( r->L2Norm() );

    // return if the residual is already low enough
    if ( current_res < terminate_tol ) {
        d_ConvergenceStatus = current_res < static_cast<T>( d_dAbsoluteTolerance ) ?
                                  SolverStatus::ConvergedOnAbsTol :
                                  SolverStatus::ConvergedOnRelTol;
        // provide history conv history etc?
        return;
    }

    z = u->clone();

    // apply the preconditioner if it exists
    if ( d_bUsesPreconditioner ) {
        d_pPreconditioner->apply( r, z );
    } else {
        z->copyVector( r );
    }

    auto rho_1 = static_cast<T>( z->dot( *r ) );
    auto rho_0 = rho_1;

    auto p = z->clone();
    auto w = r->clone();
    p->copyVector( z );

    for ( d_iNumberIterations = 0; d_iNumberIterations < d_iMaxIterations; ++d_iNumberIterations ) {

        AMP::Scalar beta{ static_cast<T>( 1.0 ) };

        p->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
        // w = Ap
        d_pOperator->apply( p, w );

        // alpha = p'Ap
        auto alpha = static_cast<T>( w->dot( *p ) );

        // sanity check, the curvature should be positive
        if ( alpha == 0.0 ) {
            // at solution
            d_ConvergenceStatus = current_res < static_cast<T>( d_dAbsoluteTolerance ) ?
                                      SolverStatus::ConvergedOnAbsTol :
                                      SolverStatus::ConvergedOnRelTol;
            break;
        } else if ( alpha < 0.0 ) {
            // set diverged reason
            d_ConvergenceStatus = SolverStatus::DivergedOther;
            break;
        }

        alpha = rho_1 / alpha;

        u->axpy( alpha, *p, *u );
        r->axpy( -alpha, *w, *r );

        // compute the current residual norm
        current_res = static_cast<T>( r->L2Norm() );
        if ( d_iDebugPrintInfoLevel > 0 ) {
            AMP::pout << "CG: iteration " << ( d_iNumberIterations + 1 ) << ", residual "
                      << current_res << ", solution norm " << u->L2Norm() << std::endl;
        }
        // check if converged
        if ( current_res < terminate_tol ) {
            d_ConvergenceStatus = current_res < static_cast<T>( d_dAbsoluteTolerance ) ?
                                      SolverStatus::ConvergedOnAbsTol :
                                      SolverStatus::ConvergedOnRelTol;
            // provide history conv history etc?
            break;
        }

        // apply the preconditioner if it exists
        if ( d_bUsesPreconditioner ) {
            d_pPreconditioner->apply( r, z );
        } else {
            z->copyVector( r );
        }

        rho_0 = rho_1;
        rho_1 = static_cast<T>( r->dot( *z ) );

        beta = rho_1 / rho_0;
        p->axpy( beta, *p, *z );
    }

    if ( d_iDebugPrintInfoLevel > 2 ) {
        AMP::pout << "L2Norm of solution: " << u->L2Norm() << std::endl;
    }

    u->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );

    if ( d_bComputeResidual ) {
        d_pOperator->residual( f, u, r );
        d_dResidualNorm = static_cast<T>( r->L2Norm() );
    } else
        d_dResidualNorm = current_res;
}

template<typename T>
void CGSolver<T>::resetOperator( std::shared_ptr<const AMP::Operator::OperatorParameters> params )
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
