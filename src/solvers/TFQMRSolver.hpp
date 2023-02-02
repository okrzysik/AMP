#include "AMP/operators/LinearOperator.h"
#include "AMP/solvers/TFQMRSolver.h"
#include "ProfilerApp.h"


#include <array>
#include <cmath>
#include <limits>

namespace AMP::Solver {

/****************************************************************
 *  Constructors                                                 *
 ****************************************************************/

template<typename T>
TFQMRSolver<T>::TFQMRSolver( std::shared_ptr<SolverStrategyParameters> parameters )
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
void TFQMRSolver<T>::initialize( std::shared_ptr<const SolverStrategyParameters> parameters )
{
    AMP_ASSERT( parameters );

    d_pPreconditioner = parameters->d_pNestedSolver;

    getFromInput( parameters->d_db );

    if ( d_pOperator ) {
        registerOperator( d_pOperator );
    }
}

// Function to get values from input
template<typename T>
void TFQMRSolver<T>::getFromInput( std::shared_ptr<const AMP::Database> db )
{
    d_bUsesPreconditioner = db->getWithDefault<bool>( "uses_preconditioner", false );

    // default is right preconditioning, options are right, left, both
    if ( d_bUsesPreconditioner ) {
        d_preconditioner_side = db->getWithDefault<std::string>( "preconditioner_side", "right" );
    }
}

/****************************************************************
 *  Solve                                                        *
 * TODO: store convergence history, iterations, convergence reason
 ****************************************************************/
template<typename T>
void TFQMRSolver<T>::apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                            std::shared_ptr<AMP::LinearAlgebra::Vector> x )
{
    PROFILE_START( "solve" );

    // Check input vector states
    AMP_ASSERT(
        ( f->getUpdateStatus() == AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED ) ||
        ( f->getUpdateStatus() == AMP::LinearAlgebra::VectorData::UpdateState::LOCAL_CHANGED ) );
    AMP_ASSERT(
        ( x->getUpdateStatus() == AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED ) ||
        ( x->getUpdateStatus() == AMP::LinearAlgebra::VectorData::UpdateState::LOCAL_CHANGED ) );

    // compute the norm of the rhs in order to compute
    // the termination criterion
    auto f_norm = static_cast<T>( f->L2Norm() );

    // if the rhs is zero we try to converge to the relative convergence
    if ( f_norm == static_cast<T>( 0.0 ) ) {
        f_norm = static_cast<T>( 1.0 );
    }

    const T terminate_tol = d_dRelativeTolerance * f_norm;

    if ( d_iDebugPrintInfoLevel > 2 ) {
        std::cout << "TFQMRSolver<T>::solve: initial L2Norm of solution vector: " << x->L2Norm()
                  << std::endl;
        std::cout << "TFQMRSolver<T>::solve: initial L2Norm of rhs vector: " << f_norm << std::endl;
    }

    if ( d_pOperator ) {
        registerOperator( d_pOperator );
    }

    // residual vector
    auto res = f->cloneVector();

    // compute the initial residual
    if ( d_bUseZeroInitialGuess ) {
        res->copyVector( f );
    } else {
        d_pOperator->residual( f, x, res );
    }

    // compute the current residual norm
    auto res_norm = static_cast<T>( res->L2Norm() );

    if ( d_iDebugPrintInfoLevel > 0 ) {
        std::cout << "TFQMR: initial residual " << res_norm << std::endl;
    }

    // return if the residual is already low enough
    if ( res_norm < terminate_tol ) {
        // provide a convergence reason
        // provide history (iterations, conv history etc)
        if ( d_iDebugPrintInfoLevel > 0 ) {
            std::cout << "TFQMRSolver<T>::solve: initial residual norm " << res_norm
                      << " is below convergence tolerance: " << terminate_tol << std::endl;
        }
        return;
    }

    // parameters in TFQMR
    T theta  = static_cast<T>( 0.0 );
    T eta    = static_cast<T>( 0.0 );
    T tau    = res_norm;
    auto rho = tau * tau;

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

        d_pPreconditioner->apply( y[0], z );

    } else {

        z = y[0];
    }

    d_pOperator->apply( z, v );

    u[0]->copyVector( v );

    int k = 0;

    bool converged = false;

    while ( k < d_iMaxIterations ) {

        ++k;
        auto sigma = static_cast<T>( res->dot( *v ) );

        // replace by soft-equal
        if ( sigma == static_cast<T>( 0.0 ) ) {
            // the method breaks down as the vectors are orthogonal to r
            AMP_ERROR( "TFQMR breakdown, sigma == 0 " );
        }

        auto alpha = rho / sigma;

        for ( int j = 0; j <= 1; ++j ) {

            if ( j == 1 ) {

                y[1]->axpy( -alpha, *v, *y[0] );

                if ( d_bUsesPreconditioner && ( d_preconditioner_side == "right" ) ) {

                    d_pPreconditioner->apply( y[1], z );

                } else {
                    z = y[1];
                }

                d_pOperator->apply( z, u[1] );
            }

            const int m = 2 * k - 1 + j;
            w->axpy( -alpha, *u[j], *w );
            d->axpy( ( theta * theta * eta / alpha ), *d, *y[j] );

            theta = static_cast<T>( w->L2Norm() ) / tau;
            const auto c =
                static_cast<T>( 1.0 ) / std::sqrt( static_cast<T>( 1.0 ) + theta * theta );
            tau = tau * theta * c;
            eta = c * c * alpha;

            // update the increment to the solution
            delta->axpy( eta, *d, *delta );

            if ( tau * ( std::sqrt( (T) ( m + 1 ) ) ) <= terminate_tol ) {

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
                d_pPreconditioner->apply( delta, z );
            } else {
                z = delta;
            }
            x->axpy( 1.0, *z, *x );
            return;
        }

        // replace by soft-equal
        if ( rho == static_cast<T>( 0.0 ) ) {
            // the method breaks down as rho==0
            AMP_ERROR( "TFQMR breakdown, rho == 0 " );
        }

        auto rho_n = static_cast<T>( res->dot( *w ) );
        auto beta  = rho_n / rho;
        rho        = rho_n;

        y[0]->axpy( beta, *y[1], *w );

        if ( d_bUsesPreconditioner && ( d_preconditioner_side == "right" ) ) {

            d_pPreconditioner->apply( y[0], z );

        } else {

            z = y[0];
        }

        d_pOperator->apply( z, u[0] );

        v->axpy( beta, *v, *u[1] );
        v->axpy( beta, *v, *u[0] );

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
template<typename T>
void TFQMRSolver<T>::registerOperator( std::shared_ptr<AMP::Operator::Operator> op )
{
    AMP_ASSERT( op );
    d_pOperator = op;
}

template<typename T>
void TFQMRSolver<T>::resetOperator(
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
