#include "AMP/solvers/TFQMRSolver.h"
#include "ProfilerApp.h"
#include "AMP/operators/LinearOperator.h"
#include "AMP/solvers/KrylovSolverParameters.h"


#include <array>
#include <cmath>
#include <limits>

namespace AMP {
namespace Solver {

/****************************************************************
 *  Constructors                                                 *
 ****************************************************************/
TFQMRSolver::TFQMRSolver() : d_restarts( 0 ) {}

TFQMRSolver::TFQMRSolver( AMP::shared_ptr<SolverStrategyParameters> parameters )
    : SolverStrategy( parameters ), d_restarts( 0 )
{
    AMP_ASSERT( parameters.get() != nullptr );

    // Initialize
    initialize( parameters );
}


/****************************************************************
 *  Destructor                                                   *
 ****************************************************************/
TFQMRSolver::~TFQMRSolver() {}

/****************************************************************
 *  Initialize                                                   *
 ****************************************************************/
void TFQMRSolver::initialize( AMP::shared_ptr<SolverStrategyParameters> const params )
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
void TFQMRSolver::getFromInput( const AMP::shared_ptr<AMP::Database> &db )
{

    d_dRelativeTolerance = db->getDoubleWithDefault( "relative_tolerance", 1.0e-9 );
    d_iMaxIterations     = db->getDoubleWithDefault( "max_iterations", 1000 );

    d_bUsesPreconditioner = db->getBoolWithDefault( "use_preconditioner", false );

    // default is right preconditioning, options are right, left, both
    if ( d_bUsesPreconditioner ) {
        d_preconditioner_side = db->getStringWithDefault( "preconditioner_side", "right" );
    }
}

/****************************************************************
 *  Solve                                                        *
 * TODO: store convergence history, iterations, convergence reason
 ****************************************************************/
void TFQMRSolver::solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                         AMP::shared_ptr<AMP::LinearAlgebra::Vector> x )
{
    PROFILE_START( "solve" );

    // Check input vector states
    AMP_ASSERT(
        ( f->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( f->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );
    AMP_ASSERT(
        ( x->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( x->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );

    // compute the norm of the rhs in order to compute
    // the termination criterion
    double f_norm = f->L2Norm();

    // if the rhs is zero we try to converge to the relative convergence
    if ( f_norm == 0.0 ) {
        f_norm = 1.0;
    }

    const double terminate_tol = d_dRelativeTolerance * f_norm;

    if ( d_iDebugPrintInfoLevel > 2 ) {
        std::cout << "TFQMRSolver::solve: initial L2Norm of solution vector: " << x->L2Norm()
                  << std::endl;
        std::cout << "TFQMRSolver::solve: initial L2Norm of rhs vector: " << f_norm << std::endl;
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
        d_pOperator->residual( f, x, res );
    }

    // compute the current residual norm
    double res_norm = res->L2Norm();

    if ( d_iDebugPrintInfoLevel > 0 ) {
        std::cout << "TFQMR: initial residual " << res_norm << std::endl;
    }

    // return if the residual is already low enough
    if ( res_norm < terminate_tol ) {
        // provide a convergence reason
        // provide history (iterations, conv history etc)
        if ( d_iDebugPrintInfoLevel > 0 ) {
            std::cout << "TFQMRSolver::solve: initial residual norm " << res_norm
                      << " is below convergence tolerance: " << terminate_tol << std::endl;
        }
        return;
    }

    // parameters in TFQMR
    double theta = 0.0;
    double eta   = 0.0;
    double tau   = res_norm;
    double rho   = tau * tau;

    std::array<AMP::LinearAlgebra::Vector::shared_ptr, 2> u;
    u[0] = f->cloneVector();
    u[1] = f->cloneVector();
    u[0]->zero();
    u[1]->zero();

    std::array<AMP::LinearAlgebra::Vector::shared_ptr, 2> y;
    y[0] = f->cloneVector();
    y[1] = f->cloneVector();
    y[0]->zero();
    y[1]->zero();

    // z is allocated only if the preconditioner is used
    AMP::LinearAlgebra::Vector::shared_ptr z;
    if ( d_bUsesPreconditioner ) {
        z = f->cloneVector();
        z->zero();
    }

    auto delta = f->cloneVector();
    delta->zero();

    auto w = res->cloneVector();
    w->copyVector( res );

    y[0]->copyVector( res );

    auto d = res->cloneVector();
    d->zero();

    auto v = res->cloneVector();

    if ( d_bUsesPreconditioner && ( d_preconditioner_side == "right" ) ) {

        d_pPreconditioner->solve( y[0], z );

    } else {

        z = y[0];
    }

    d_pOperator->apply( z, v );

    u[0]->copyVector( v );

    int k = 0;

    bool converged = false;

    while ( k < d_iMaxIterations ) {

        ++k;
        auto sigma = res->dot( v );

        // replace by soft-equal
        if ( sigma == 0.0 ) {
            // the method breaks down as the vectors are orthogonal to r
            AMP_ERROR( "TFQMR breakdown, sigma == 0 " );
        }

        auto alpha = rho / sigma;

        for ( int j = 0; j <= 1; ++j ) {

            if ( j == 1 ) {

                y[1]->axpy( -alpha, v, y[0] );

                if ( d_bUsesPreconditioner && ( d_preconditioner_side == "right" ) ) {

                    d_pPreconditioner->solve( y[1], z );

                } else {
                    z = y[1];
                }

                d_pOperator->apply( z, u[1] );
            }

            const int m = 2 * k - 1 + j;
            w->axpy( -alpha, u[j], w );
            d->axpy( ( theta * theta * eta / alpha ), d, y[j] );

            theta        = w->L2Norm() / tau;
            const auto c = 1.0 / std::sqrt( 1 + theta * theta );
            tau          = tau * theta * c;
            eta          = c * c * alpha;

            // update the increment to the solution
            delta->axpy( eta, d, delta );

            if ( tau * ( std::sqrt( (double) ( m + 1 ) ) ) <= terminate_tol ) {

                if ( d_iDebugPrintInfoLevel > 0 ) {
                    std::cout << "TFQMR: iteration " << ( k + 1 ) << ", residual " << tau
                              << std::endl;
                }

                converged = true;
                break;
            }
        }

        if ( converged ) {

            // unwind the preconditioner if necessary
            if ( d_bUsesPreconditioner && ( d_preconditioner_side == "right" ) ) {
                d_pPreconditioner->solve( delta, z );
            } else {
                z = delta;
            }
            x->axpy( 1.0, z, x );
            return;
        }

        // replace by soft-equal
        if ( rho == 0.0 ) {
            // the method breaks down as rho==0
            AMP_ERROR( "TFQMR breakdown, rho == 0 " );
        }

        const auto rho_n = res->dot( w );
        const auto beta  = rho_n / rho;
        rho              = rho_n;

        y[0]->axpy( beta, y[1], w );

        if ( d_bUsesPreconditioner && ( d_preconditioner_side == "right" ) ) {

            d_pPreconditioner->solve( y[0], z );

        } else {

            z = y[0];
        }

        d_pOperator->apply( z, u[0] );

        v->axpy( beta, v, u[1] );
        v->axpy( beta, v, u[0] );

        if ( d_iDebugPrintInfoLevel > 0 ) {
            std::cout << "TFQMR: iteration " << ( k + 1 ) << ", residual " << tau << std::endl;
        }
    }

    if ( d_iDebugPrintInfoLevel > 2 ) {
        std::cout << "L2Norm of solution: " << x->L2Norm() << std::endl;
    }

    PROFILE_STOP( "solve" );
}

/****************************************************************
 *  Function to set the register the operator                    *
 ****************************************************************/
void TFQMRSolver::registerOperator( const AMP::shared_ptr<AMP::Operator::Operator> op )
{
    AMP_ASSERT( op.get() != nullptr );

    d_pOperator = op;

    AMP::shared_ptr<AMP::Operator::LinearOperator> linearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( op );
    AMP_ASSERT( linearOperator.get() != nullptr );
}
void TFQMRSolver::resetOperator( const AMP::shared_ptr<AMP::Operator::OperatorParameters> params )
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
