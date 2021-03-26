#include "AMP/solvers/GMRESSolver.h"
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
GMRESSolver::GMRESSolver() : d_restarts( 0 ) { NULL_USE( d_restarts ); }

GMRESSolver::GMRESSolver( std::shared_ptr<SolverStrategyParameters> parameters )
    : SolverStrategy( parameters ), d_restarts( 0 )
{
    AMP_ASSERT( parameters.get() != nullptr );

    // Initialize
    initialize( parameters );
}


/****************************************************************
 *  Destructor                                                   *
 ****************************************************************/
GMRESSolver::~GMRESSolver() {}

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

    // maximum dimension to allocate storage for
    const int max_dim = std::min( d_iMaxKrylovDimension, d_iMaxIterations );
    d_dHessenberg.resize( max_dim + 1, max_dim + 1 );
    d_dHessenberg.fill( 0.0 );

    d_dcos.resize( max_dim + 1, 0.0 );
    d_dsin.resize( max_dim + 1, 0.0 );
    d_dw.resize( max_dim + 1, 0.0 );
    d_dy.resize( max_dim, 0.0 );

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

    d_sOrthogonalizationMethod = db->getWithDefault<std::string>( "ortho_method", "MGS" );

    // default is right preconditioning, options are right, left, both
    if ( d_bUsesPreconditioner ) {
        d_preconditioner_side = db->getWithDefault<std::string>( "preconditioner_side", "right" );
    }

    d_bRestart            = db->getWithDefault<bool>( "gmres_restart", false );
    d_bUsesPreconditioner = db->getWithDefault<bool>( "uses_preconditioner", false );
    d_bRestart            = db->getWithDefault<bool>( "gmres_restart", false );
}

/****************************************************************
 *  Solve                                                        *
 * TODO: store convergence history, iterations, convergence reason
 ****************************************************************/
void GMRESSolver::apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
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
    // NOTE:: update this test for a better 'almost equal'
    if ( f_norm < std::numeric_limits<double>::epsilon() ) {
        f_norm = 1.0;
    }

    const double terminate_tol = d_dRelativeTolerance * f_norm;

    if ( d_iDebugPrintInfoLevel > 2 ) {
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
        u->setToScalar( 0.0 );
    } else {
        d_pOperator->residual( f, u, res );
    }

    d_nr = -1;

    // compute the current residual norm
    const double beta = static_cast<double>( res->L2Norm() );

    if ( d_iDebugPrintInfoLevel > 0 ) {
        std::cout << "GMRES: initial residual " << beta << std::endl;
    }

    // return if the residual is already low enough
    if ( beta < terminate_tol ) {
        if ( d_iDebugPrintInfoLevel > 0 ) {
            std::cout << "GMRESSolver::solve: initial residual norm " << beta
                      << " is below convergence tolerance: " << terminate_tol << std::endl;
        }

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

    // z is only used if there is preconditioning
    AMP::LinearAlgebra::Vector::shared_ptr z;

    for ( int k = 0; ( k < d_iMaxIterations ) && ( v_norm > terminate_tol ); ++k ) {

        // clone off of the rhs to create a new basis vector
        AMP::LinearAlgebra::Vector::shared_ptr v = f->cloneVector();
        if ( d_bUsesPreconditioner && ( d_preconditioner_side == "right" ) ) {
            z = f->cloneVector();
            d_pPreconditioner->apply( d_vBasis[k], z );
        } else {
            z = d_vBasis[k];
        }

        // construct the Krylov vector
        d_pOperator->apply( z, v );

        // orthogonalize to previous vectors and
        // add new column to Hessenberg matrix
        orthogonalize( v );

        v_norm = d_dHessenberg( k + 1, k );
        // replace the conditional by a soft equality
        // check for happy breakdown
        if ( v_norm != 0.0 ) {
            v->scale( 1.0 / v_norm );
            v->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
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
            d_dHessenberg( k + 1, k ) = 0.0; // explicitly set subdiag to zero to prevent round-off

            // explicitly apply the newly computed
            // Givens rotations to the rhs vector
            auto x = d_dw[k];
            auto y = d_dw[k + 1];
            auto c = d_dcos[k];
            auto s = d_dsin[k];
#if 0
            d_dw[k]     = c * x - s * y;
            d_dw[k + 1] = s * x + c * y;
#else
            d_dw[k]     = c * x + s * y;
            d_dw[k + 1] = -s * x + c * y;
#endif
        }

        v_norm = std::fabs( d_dw[k + 1] );

        if ( d_iDebugPrintInfoLevel > 0 ) {
            std::cout << "GMRES: iteration " << ( k + 1 ) << ", residual " << v_norm << std::endl;
        }

        ++d_nr; // update the dimension of the upper triangular system to solve
    }

    // compute y, the solution to the least squares minimization problem
    backwardSolve();

    // update the current approximation with the correction
    if ( d_bUsesPreconditioner && ( d_preconditioner_side == "right" ) ) {

        z->setToScalar( 0.0 );

        for ( int i = 0; i <= d_nr; ++i ) {
            z->axpy( d_dy[i], *d_vBasis[i], *z );
        }

        AMP::LinearAlgebra::Vector::shared_ptr v = f->cloneVector();
        d_pPreconditioner->apply( z, v );
        u->axpy( 1.0, *v, *u );

    } else {
        for ( int i = 0; i <= d_nr; ++i ) {
            u->axpy( d_dy[i], *d_vBasis[i], *u );
        }
    }
    u->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    if ( d_iDebugPrintInfoLevel > 2 ) {
        d_pOperator->residual( f, u, res );
        std::cout << "GMRES: Final residual: " << res->L2Norm() << std::endl;
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

            const double h_jk = static_cast<double>( v->dot( *d_vBasis[j] ) );
            v->axpy( -h_jk, *d_vBasis[j], *v );
            d_dHessenberg( j, k - 1 ) = h_jk;
        }
    } else {

        AMP_ERROR( "Unknown orthogonalization method in GMRES" );
    }

    v->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    // h_{k+1, k}
    const auto v_norm         = static_cast<double>( v->L2Norm() );
    d_dHessenberg( k, k - 1 ) = v_norm; // adjusting for zero starting index
}

void GMRESSolver::applyGivensRotation( const int i, const int k )
{
    // updates column k of the Hessenberg matrix by applying the i-th Givens rotations

    auto x = d_dHessenberg( i, k );
    auto y = d_dHessenberg( i + 1, k );
    auto c = d_dcos[i];
    auto s = d_dsin[i];

#if 0
    d_dHessenberg( i, k )     = c * x - s * y;
    d_dHessenberg( i + 1, k ) = s * x + c * y;
#else
    d_dHessenberg( i, k ) = c * x + s * y;
    d_dHessenberg( i + 1, k ) = -s * x + c * y;
#endif
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

    decltype( f ) c, s;
    if ( g == 0.0 ) {
        c = 1.0;
        s = 0.0;
    } else if ( f == 0.0 ) {
        c = 0.0;
        s = ( g < 0.0 ) ? -1.0 : 1.0;
    } else {
        decltype( f ) r;
        r = std::sqrt( f * f + g * g );
        r = 1.0 / r;
        c = std::fabs( f ) * r;
        s = std::copysign( g * r, f );
    }

    d_dcos[k] = c;
    d_dsin[k] = s;
}

void GMRESSolver::backwardSolve( void )
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

    auto linearOperator = std::dynamic_pointer_cast<AMP::Operator::LinearOperator>( op );
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
