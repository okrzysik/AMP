#include "AMP/solvers/GMRESSolver.h"
#include "AMP/operators/LinearOperator.h"
#include "ProfilerApp.h"


#include <cmath>
#include <limits>

namespace AMP {
namespace Solver {

/****************************************************************
 *  Constructors                                                 *
 ****************************************************************/
GMRESSolver::GMRESSolver() : d_restarts( 0 ) { NULL_USE( d_restarts ); }

GMRESSolver::GMRESSolver( std::shared_ptr<KrylovSolverParameters> parameters )
    : SolverStrategy( parameters ), d_restarts( 0 )
{
    AMP_ASSERT( parameters.get() != nullptr );

    // Initialize
    initialize( parameters );
}


/****************************************************************
 *  Destructor                                                   *
 ****************************************************************/
GMRESSolver::~GMRESSolver() = default;

/****************************************************************
 *  Initialize                                                   *
 ****************************************************************/
void GMRESSolver::initialize( std::shared_ptr<SolverStrategyParameters> const params )
{
    auto parameters = std::dynamic_pointer_cast<KrylovSolverParameters>( params );
    AMP_ASSERT( parameters.get() != nullptr );
    d_comm = parameters->d_comm;
    AMP_ASSERT( !d_comm.isNull() );

    getFromInput( parameters->d_db );

    d_dHessenberg.resize( d_iMaxKrylovDimension + 1, d_iMaxKrylovDimension + 1 );
    d_dHessenberg.fill( 0.0 );

    d_dcos.resize( d_iMaxKrylovDimension + 1, 0.0 );
    d_dsin.resize( d_iMaxKrylovDimension + 1, 0.0 );
    d_dw.resize( d_iMaxKrylovDimension + 1, 0.0 );
    d_dy.resize( d_iMaxKrylovDimension, 0.0 );

    d_pPreconditioner = parameters->d_pPreconditioner;

    if ( d_pOperator.get() != nullptr ) {
        registerOperator( d_pOperator );
    }
}

// Function to get values from input
void GMRESSolver::getFromInput( std::shared_ptr<AMP::Database> db )
{
    // the max iterations could be larger than the max Krylov dimension
    // in the case of restarted GMRES so we allow specification separately
    d_iMaxKrylovDimension = db->getWithDefault<double>( "max_dimension", 100 );
    d_iMaxIterations      = db->getWithDefault<double>( "max_iterations", d_iMaxKrylovDimension );

    d_dRelativeTolerance = db->getWithDefault<double>( "relative_tolerance", 1.0e-9 );

    d_sOrthogonalizationMethod = db->getWithDefault<std::string>( "ortho_method", "MGS" );

    d_bUsesPreconditioner = db->getWithDefault( "uses_preconditioner", false );
    d_bRestart            = db->getWithDefault( "gmres_restart", false );
}

/****************************************************************
 *  Solve                                                        *
 * TODO: store convergence history, iterations, convergence reason
 ****************************************************************/
void GMRESSolver::solve( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
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
    double f_norm = f->L2Norm();

    // if the rhs is zero we try to converge to the relative convergence
    if ( f_norm == 0.0 ) {
        f_norm = 1.0;
    }

    const double terminate_tol = d_dRelativeTolerance * f_norm;

    if ( d_iDebugPrintInfoLevel > 1 ) {
        std::cout << "GMRESSolver::solve: initial L2Norm of solution vector: " << u->L2Norm()
                  << std::endl;
        std::cout << "GMRESSolver::solve: initial L2Norm of rhs vector: " << f_norm << std::endl;
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
    const double beta = res->L2Norm();

    // return if the residual is already low enough
    if ( beta < terminate_tol ) {
        // provide a convergence reason
        // provide history (iterations, conv history etc)
        return;
    }

    // normalize the first basis vector
    res->scale( 1.0 / beta );

    // push the residual as the first basis vector
    d_vBasis.push_back( res );

    // 'w*e_1' is the rhs for the least squares minimization problem
    d_dw[0] = beta;

    auto v_norm = beta;

    for ( int k = 0; ( k < d_iMaxIterations ) && ( v_norm < terminate_tol ); ++k ) {

        // clone off of the rhs to create a new basis vector
        AMP::LinearAlgebra::Vector::shared_ptr v = f->cloneVector();

        // construct the Krylov vector
        d_pOperator->apply( d_vBasis[k], v );

        // orthogonalize to previous vectors and
        // add new column to Hessenberg matrix
        orthogonalize( v );

        v_norm = d_dHessenberg( k + 1, k );
        // replace the conditional by a soft equality
        // check for happy breakdown
        if ( v_norm != 0.0 ) {
            v->scale( 1.0 / v_norm );
        }

        // update basis with new orthonormal vector
        d_vBasis.push_back( v );

        // apply all previous Givens rotations to
        // the k-th column of the Hessenberg matrix
        for ( int i = 0; i < k; ++i ) {
            applyGivensRotation( i, k );
        }

        if ( v_norm != 0.0 ) {
            // compute and store the Givens rotation that zeroes out
            // the subdiagonal for the current column
            computeGivensRotation( k );
            // zero out the subdiagonal
            applyGivensRotation( k, k );

            // explicitly apply the newly computed
            // Givens rotations to the rhs vector
            auto x      = d_dw[k];
            auto y      = d_dw[k];
            auto c      = d_dcos[k];
            auto s      = d_dsin[k];
            d_dw[k]     = c * x - s * y;
            d_dw[k + 1] = s * x + c * y;
        }

        // check on whether the comparison should be to terminate_tol
        if ( std::fabs( d_dw[k + 1] ) < terminate_tol ) {
            d_nr = k; // set the dimension of the upper triangular system to solve
            break;
        }
    }

    // compute y, the solution to the least squares minimization problem
    backwardSolve();

    // update the current approximation with the correction
    for ( int i = 0; i <= d_nr; ++i ) {
        u->axpy( d_dy[i], d_vBasis[i], u );
    }

    // compute a final residual ??
    d_pOperator->residual( f, u, res );

    if ( d_iDebugPrintInfoLevel > 2 ) {
        std::cout << "L2Norm of solution: " << u->L2Norm() << std::endl;
    }

    PROFILE_STOP( "solve" );
}

void GMRESSolver::orthogonalize( std::shared_ptr<AMP::LinearAlgebra::Vector> v )
{
    const int k = d_vBasis.size();

    if ( d_sOrthogonalizationMethod == "CGS" ) {

        AMP_ERROR( "Classical Gram-Schmidt not implemented as yet" );
    } else if ( d_sOrthogonalizationMethod == "MGS" ) {

        for ( int j = 0; j < k; ++j ) {

            const double h_jk = v->dot( d_vBasis[j] );
            v->axpy( -h_jk, d_vBasis[j], v );
            d_dHessenberg( j, k ) = h_jk;
        }
    } else {

        AMP_ERROR( "Unknown orthogonalization method in GMRES" );
    }

    // h_{k+1, k}
    auto v_norm               = v->L2Norm();
    d_dHessenberg( k, k - 1 ) = v_norm; // adjusting for zero starting index
}

void GMRESSolver::applyGivensRotation( const int i, const int k )
{
    // updates column k of the Hessenberg matrix by applying the i-th Givens rotations

    auto x = d_dHessenberg( i, k );
    auto y = d_dHessenberg( i + 1, k );
    auto c = d_dcos[i];
    auto s = d_dsin[i];

    d_dHessenberg( i, k )     = c * x - s * y;
    d_dHessenberg( i + 1, k ) = s * x + c * y;
}

void GMRESSolver::computeGivensRotation( const int k )
{
    // computes the Givens rotation required to zero out
    // the subdiagonal on column k of the Hessenberg matrix

    // The implementation here follows Algorithm 1 in
    // "On Computing Givens rotations reliably and efficiently"
    // by D. Bindel, J. Demmel, W. Kahan, O. Marques
    // UT-CS-00-449, October 2000.

    auto f = d_dHessenberg( k, k );
    auto g = d_dHessenberg( k + 1, k );

    auto c = f;
    auto s = c;
    auto r = c;
    NULL_USE( r );

    if ( g == 0.0 ) {

        c = 1.0;
        s = 0.0;
    } else if ( f == 0.0 ) {

        c = 0.0;
        s = ( g < 0.0 ) ? -1.0 : 1.0;
    } else {

        r = std::sqrt( f * f + g * g );
        c = std::fabs( f ) / r;
        s = std::copysign( g / r, f );
    }

    d_dcos[k] = c;
    d_dsin[k] = s;
}

void GMRESSolver::backwardSolve()
{
    // lower corner
    d_dy[d_nr] = d_dw[d_nr] / d_dHessenberg( d_nr, d_nr );

    // backwards solve
    for ( int k = d_nr - 1; k >= 0; --k ) {

        d_dy[k] = d_dw[k];

        for ( int i = k + 1; i <= d_nr; ++i ) {
            d_dy[k] -= d_dHessenberg( k, i ) * d_dy[i];
        }

        d_dy[k] = d_dy[k] / d_dHessenberg( k, k );
    }
}

/****************************************************************
 *  Function to set the register the operator                    *
 ****************************************************************/
void GMRESSolver::registerOperator( const std::shared_ptr<AMP::Operator::Operator> op )
{
    AMP_ASSERT( op.get() != nullptr );

    d_pOperator = op;

    std::shared_ptr<AMP::Operator::LinearOperator> linearOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearOperator>( op );
    AMP_ASSERT( linearOperator.get() != nullptr );
}
void GMRESSolver::resetOperator( const std::shared_ptr<AMP::Operator::OperatorParameters> params )
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
