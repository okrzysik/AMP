#include "AMP/solvers/CGSolver.h"
#include "AMP/operators/LinearOperator.h"
#include "AMP/solvers/KrylovSolverParameters.h"
#include "ProfilerApp.h"


namespace AMP {
namespace Solver {

/****************************************************************
 *  Constructors                                                 *
 ****************************************************************/
CGSolver::CGSolver() {}

CGSolver::CGSolver( AMP::shared_ptr<SolverStrategyParameters> parameters )
    : SolverStrategy( parameters )
{
    AMP_ASSERT( parameters.get() != nullptr );

    // Initialize
    initialize( parameters );
}


/****************************************************************
 *  Destructor                                                   *
 ****************************************************************/
CGSolver::~CGSolver() {}

/****************************************************************
 *  Initialize                                                   *
 ****************************************************************/
void CGSolver::initialize( AMP::shared_ptr<SolverStrategyParameters> const params )
{
    auto parameters = AMP::dynamic_pointer_cast<KrylovSolverParameters>( params );
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
void CGSolver::getFromInput( const AMP::shared_ptr<AMP::Database> &db )
{

    d_dRelativeTolerance   = db->getDoubleWithDefault( "relative_tolerance", 1.0e-9 );
    d_dAbsoluteTolerance   = db->getDoubleWithDefault( "absolute_tolerance", 1.0e-14 );
    d_dDivergenceTolerance = db->getDoubleWithDefault( "divergence_tolerance", 1.0e+03 );
    d_iMaxIterations       = db->getDoubleWithDefault( "max_iterations", 1000 );

    d_bUsesPreconditioner = db->getBoolWithDefault( "uses_preconditioner", false );
}

/****************************************************************
 *  Solve                                                        *
 * TODO: store convergence history, iterations, convergence reason
 ****************************************************************/
void CGSolver::solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                      AMP::shared_ptr<AMP::LinearAlgebra::Vector> u )
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
    const double f_norm = f->L2Norm();

    // enhance with convergence reason, number of iterations etc
    if ( f_norm == 0.0 )
        return;

    const double terminate_tol = d_dRelativeTolerance * f_norm;

    if ( d_iDebugPrintInfoLevel > 1 ) {
        std::cout << "CGSolver::solve: initial L2Norm of solution vector: " << u->L2Norm()
                  << std::endl;
        std::cout << "CGSolver::solve: initial L2Norm of rhs vector: " << f_norm << std::endl;
    }

    if ( d_pOperator.get() != nullptr ) {
        registerOperator( d_pOperator );
    }

    // z will store r when a preconditioner is not present
    // and will store the result of a preconditioner solve
    // when a preconditioner is present
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> z;

    // residual vector
    AMP::LinearAlgebra::Vector::shared_ptr r = f->cloneVector();

    // compute the initial residual
    if ( d_bUseZeroInitialGuess ) {
        r->copyVector( f );
    } else {
        d_pOperator->residual( f, u, r );
    }
    // compute the current residual norm
    double current_res = r->L2Norm();

    // return if the residual is already low enough
    if ( current_res < terminate_tol ) {
        // provide a convergence reason
        // provide history (iterations, conv history etc)
        return;
    }

    z = u->cloneVector();

    // apply the preconditioner if it exists
    if ( d_bUsesPreconditioner ) {
        d_pPreconditioner->solve( r, z );
    } else {
        z->copyVector( r );
    }

    std::vector<double> rho( 2, 0.0 );
    rho[1] = z->dot( r );
    rho[0] = rho[1];

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> p = z->cloneVector();
    p->copyVector( z );

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> w;

    for ( auto iter = 0; iter < d_iMaxIterations; ++iter ) {

        double beta = 1.0;

        // w = Ap
        d_pOperator->apply( p, w );

        // alpha = p'Ap
        double alpha = w->dot( p );

        // sanity check, the curvature should be positive
        if ( alpha <= 0.0 ) {
            // set diverged reason
            AMP_ERROR( "Negative curvature encountered in CG!!" );
        }

        alpha = rho[1] / alpha;

        u->axpy( alpha, p, u );
        r->axpy( -alpha, w, r );

        // compute the current residual norm
        current_res = r->L2Norm();
        // check if converged
        if ( current_res < terminate_tol ) {
            // set a convergence reason
            break;
        }

        // apply the preconditioner if it exists
        if ( d_bUsesPreconditioner ) {
            d_pPreconditioner->solve( r, z );
        } else {
            z->copyVector( r );
        }

        rho[0] = rho[1];
        rho[1] = r->dot( z );

        beta = rho[1] / rho[0];
        p->axpy( beta, p, z );
    }

    if ( d_iDebugPrintInfoLevel > 2 ) {
        std::cout << "L2Norm of solution: " << u->L2Norm() << std::endl;
    }

    PROFILE_STOP( "solve" );
}

/****************************************************************
 *  Function to set the register the operator                    *
 ****************************************************************/
void CGSolver::registerOperator( const AMP::shared_ptr<AMP::Operator::Operator> op )
{
    AMP_ASSERT( op.get() != nullptr );

    d_pOperator = op;

    AMP::shared_ptr<AMP::Operator::LinearOperator> linearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( op );
    AMP_ASSERT( linearOperator.get() != nullptr );
}
void CGSolver::resetOperator( const AMP::shared_ptr<AMP::Operator::OperatorParameters> params )
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
