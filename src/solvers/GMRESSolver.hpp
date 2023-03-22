#include "AMP/operators/LinearOperator.h"
#include "AMP/solvers/GMRESSolver.h"
#include "ProfilerApp.h"


#include <cmath>
#include <limits>

namespace AMP::Solver {

/****************************************************************
 *  Constructors                                                 *
 ****************************************************************/
template<typename T>
GMRESSolver<T>::GMRESSolver() : d_restarts( 0 )
{
    NULL_USE( d_restarts );
}

template<typename T>
GMRESSolver<T>::GMRESSolver( std::shared_ptr<SolverStrategyParameters> parameters )
    : SolverStrategy( parameters ), d_restarts( 0 )
{
    AMP_ASSERT( parameters );

    // Initialize
    initialize( parameters );
}

/****************************************************************
 *  Initialize                                                   *
 ****************************************************************/
template<typename T>
void GMRESSolver<T>::initialize( std::shared_ptr<const SolverStrategyParameters> parameters )
{
    AMP_ASSERT( parameters );

    getFromInput( parameters->d_db );

    // maximum dimension to allocate storage for
    const int max_dim = std::min( d_iMaxKrylovDimension, d_iMaxIterations );
    d_dHessenberg.resize( max_dim + 1, max_dim + 1 );
    d_dHessenberg.fill( 0.0 );

    d_dcos.resize( max_dim + 1, 0.0 );
    d_dsin.resize( max_dim + 1, 0.0 );
    d_dw.resize( max_dim + 1, 0.0 );
    d_dy.resize( max_dim, 0.0 );

    d_pPreconditioner = parameters->d_pNestedSolver;

    if ( d_pOperator ) {
        registerOperator( d_pOperator );
    }
}

// Function to get values from input
template<typename T>
void GMRESSolver<T>::getFromInput( std::shared_ptr<AMP::Database> db )
{
    // the max iterations could be larger than the max Krylov dimension
    // in the case of restarted GMRES so we allow specification separately
    d_iMaxKrylovDimension      = db->getWithDefault<int>( "max_dimension", 100 );
    d_sOrthogonalizationMethod = db->getWithDefault<std::string>( "ortho_method", "MGS" );

    // default is right preconditioning, options are right, left, both
    d_bUsesPreconditioner = db->getWithDefault<bool>( "uses_preconditioner", false );
    if ( d_bUsesPreconditioner ) {
        d_preconditioner_side = db->getWithDefault<std::string>( "preconditioner_side", "right" );
    }

    d_bRestart = db->getWithDefault<bool>( "gmres_restart", false );

    d_bFlexibleGMRES = db->getWithDefault<bool>( "flexible_gmres", false );

    if ( !d_bUsesPreconditioner )
        d_bFlexibleGMRES = false;

    if ( d_bUsesPreconditioner && d_preconditioner_side == "left" && d_bFlexibleGMRES )
        AMP_ERROR( "Flexible GMRES needs right preconditioning" );
}

/****************************************************************
 *  Solve                                                        *
 * TODO: store convergence history, iterations, convergence reason
 ****************************************************************/
template<typename T>
void GMRESSolver<T>::apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
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
    auto f_norm = static_cast<T>( f->L2Norm() );

    // if the rhs is zero we try to converge to the relative convergence
    // NOTE:: update this test for a better 'almost equal'
    if ( f_norm < std::numeric_limits<T>::epsilon() ) {
        f_norm = static_cast<T>( 1.0 );
    }

    const T terminate_tol = d_dRelativeTolerance * f_norm;

    if ( d_iDebugPrintInfoLevel > 2 ) {
        std::cout << "GMRESSolver<T>::solve: initial L2Norm of solution vector: " << u->L2Norm()
                  << std::endl;
        std::cout << "GMRESSolver<T>::solve: initial L2Norm of rhs vector: " << f_norm << std::endl;
    }

    if ( d_pOperator ) {
        registerOperator( d_pOperator );
    }

    // residual vector
    AMP::LinearAlgebra::Vector::shared_ptr res = f->clone();

    // z is only used if there is preconditioning
    AMP::LinearAlgebra::Vector::shared_ptr z;
    // z1 is only used if there is right preconditioning
    AMP::LinearAlgebra::Vector::shared_ptr z1;

    if ( d_bUsesPreconditioner ) {
        z = f->clone();
        if ( d_preconditioner_side == "right" ) {
            z1 = f->clone();
        }
    }

    // compute the initial residual
    computeInitialResidual( d_bUseZeroInitialGuess, f, u, z, res );

    // compute the current residual norm
    const auto beta = static_cast<T>( res->L2Norm() );

    if ( d_iDebugPrintInfoLevel > 0 ) {
        std::cout << "GMRES: initial residual " << beta << std::endl;
    }

    // return if the residual is already low enough
    if ( beta < terminate_tol ) {
        if ( d_iDebugPrintInfoLevel > 0 ) {
            std::cout << "GMRESSolver<T>::solve: initial residual norm " << beta
                      << " is below convergence tolerance: " << terminate_tol << std::endl;
        }

        // provide history (iterations, conv history etc)
        return;
    }

    // normalize the first basis vector
    res->scale( static_cast<T>( 1.0 ) / beta );

    // push the residual as the first basis vector
    d_vBasis.push_back( res );

    // 'w*e_1' is the rhs for the least squares minimization problem
    d_dw[0] = beta;

    auto v_norm = beta;

    int k = 0;
    for ( int iter = 0; ( iter < d_iMaxIterations ) && ( v_norm > terminate_tol ); ++iter ) {

        AMP::LinearAlgebra::Vector::shared_ptr v;
        AMP::LinearAlgebra::Vector::shared_ptr zb;

        if ( k + 1 < static_cast<int>( d_vBasis.size() ) ) {
            // reuse basis vectors for restarts in order to avoid a clone
            v = d_vBasis[k + 1];
            if ( d_bFlexibleGMRES )
                // z_Basis is increased in size one behind wrt d_vBasis
                zb = d_zBasis[k];
        } else {
            // clone off of the rhs to create a new basis vector
            v = f->clone();
            d_vBasis.push_back( v );
            if ( d_bFlexibleGMRES )
                zb = f->cloneVector();
        }

        if ( d_bUsesPreconditioner && ( d_preconditioner_side == "left" ) ) {
            d_pOperator->apply( d_vBasis[k], z );
            // construct the Krylov vector
            d_pPreconditioner->apply( z, v );
        } else {
            if ( d_bUsesPreconditioner && ( d_preconditioner_side == "right" ) ) {
                if ( !d_bFlexibleGMRES ) {
                    d_pPreconditioner->apply( d_vBasis[k], z );
                    d_pOperator->apply( z, v );
                } else {
                    d_pPreconditioner->apply( d_vBasis[k], zb );
                    d_zBasis.push_back( zb );
                    d_pOperator->apply( zb, v );
                }
            } else {
                z = d_vBasis[k];
                d_pOperator->apply( z, v );
            }
        }

        // orthogonalize to previous vectors and
        // add new column to Hessenberg matrix
        orthogonalize( k + 1, v );

        v_norm = d_dHessenberg( k + 1, k );
        // replace the conditional by a soft equality
        // check for happy breakdown
        if ( v_norm != static_cast<T>( 0.0 ) ) {
            v->scale( static_cast<T>( 1.0 ) / v_norm );
            v->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
        }

        // apply all previous Givens rotations to
        // the k-th column of the Hessenberg matrix
        for ( int i = 0; i < k; ++i ) {
            applyGivensRotation( i, k );
        }

        if ( v_norm != static_cast<T>( 0.0 ) ) {
            // compute and store the Givens rotation that zeroes out
            // the subdiagonal for the current column
            computeGivensRotation( k );
            // zero out the subdiagonal
            applyGivensRotation( k, k );
            // explicitly set subdiag to zero to prevent round-off
            d_dHessenberg( k + 1, k ) = static_cast<T>( 0.0 );

            // explicitly apply the newly computed
            // Givens rotations to the rhs vector
            auto x = d_dw[k];
            auto c = d_dcos[k];
            auto s = d_dsin[k];

            d_dw[k]     = c * x;
            d_dw[k + 1] = -s * x;
        }

        // this is the norm of the residual thanks to the Givens rotations
        v_norm = std::fabs( d_dw[k + 1] );

        if ( d_iDebugPrintInfoLevel > 0 ) {
            std::cout << "GMRES: iteration " << ( iter + 1 ) << ", residual " << v_norm
                      << std::endl;
        }

        ++k;
        if ( ( k == d_iMaxKrylovDimension ) && ( iter != d_iMaxIterations - 1 ) ) {
            if ( d_bRestart ) {
                // with restarts, you start over with the last solution as new initial guess to
                // compute the residal: r^new = Ax^old - b
                if ( d_iDebugPrintInfoLevel > 2 ) {
                    std::cout << "GMRES: restarting" << std::endl;
                }
                // note: backwardSolve and addCorrection only go up to k - 1 because k has already
                // been increased with ++k (so they really go up to k)
                // compute y, the solution to the least squares minimization problem
                backwardSolve( k - 1 );
                // update the current approximation with the correction
                addCorrection( k - 1, z, z1, u );
                computeInitialResidual( false, f, u, z, d_vBasis[0] );
                v_norm = static_cast<T>( d_vBasis[0]->L2Norm() );
                d_vBasis[0]->scale( static_cast<T>( 1.0 ) / v_norm );
                d_dw[0] = v_norm;
                ++d_restarts;
                k = 0;
            } else
                break;
        }
    }

    if ( k > 0 ) {
        // compute y, the solution to the least squares minimization problem
        backwardSolve( k - 1 );

        // update the current approximation with the correction
        addCorrection( k - 1, z, z1, u );
    }
    if ( d_iDebugPrintInfoLevel > 2 ) {
        d_pOperator->residual( f, u, res );
        std::cout << "GMRES: Final residual: " << res->L2Norm() << std::endl;
        std::cout << "L2Norm of solution: " << u->L2Norm() << std::endl;
    }

    PROFILE_STOP( "solve" );
}

template<typename T>
void GMRESSolver<T>::orthogonalize( const int k, std::shared_ptr<AMP::LinearAlgebra::Vector> v )
{
    if ( d_sOrthogonalizationMethod == "CGS" ) {

        AMP_ERROR( "Classical Gram-Schmidt not implemented as yet" );
    } else if ( d_sOrthogonalizationMethod == "MGS" ) {

        for ( int j = 0; j < k; ++j ) {

            const auto h_jk = static_cast<T>( v->dot( *d_vBasis[j] ) );
            v->axpy( -h_jk, *d_vBasis[j], *v );
            d_dHessenberg( j, k - 1 ) = h_jk;
        }
    } else {

        AMP_ERROR( "Unknown orthogonalization method in GMRES" );
    }

    v->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    // h_{k+1, k}
    const auto v_norm         = static_cast<T>( v->L2Norm() );
    d_dHessenberg( k, k - 1 ) = v_norm; // adjusting for zero starting index
}

template<typename T>
void GMRESSolver<T>::applyGivensRotation( const int i, const int k )
{
    // updates column k of the Hessenberg matrix by applying the i-th Givens rotations

    auto x = d_dHessenberg( i, k );
    auto y = d_dHessenberg( i + 1, k );
    auto c = d_dcos[i];
    auto s = d_dsin[i];

    d_dHessenberg( i, k )     = c * x + s * y;
    d_dHessenberg( i + 1, k ) = -s * x + c * y;
}

template<typename T>
void GMRESSolver<T>::computeGivensRotation( const int k )
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
    if ( g == static_cast<T>( 0.0 ) ) {
        c = static_cast<T>( 1.0 );
        s = static_cast<T>( 0.0 );
    } else if ( f == static_cast<T>( 0.0 ) ) {
        c = static_cast<T>( 0.0 );
        s = ( g < static_cast<T>( 0.0 ) ) ? static_cast<T>( -1.0 ) : static_cast<T>( 1.0 );
    } else {
        decltype( f ) r;
        r = std::sqrt( f * f + g * g );
        r = static_cast<T>( 1.0 ) / r;
        c = std::fabs( f ) * r;
        s = std::copysign( g * r, f );
    }

    d_dcos[k] = c;
    d_dsin[k] = s;
}

template<typename T>
void GMRESSolver<T>::backwardSolve( const int nr )
{
    // lower corner
    d_dy[nr] = d_dw[nr] / d_dHessenberg( nr, nr );

    // backwards solve
    for ( int k = nr - 1; k >= 0; --k ) {

        d_dy[k] = d_dw[k];

        for ( int i = k + 1; i <= nr; ++i ) {
            d_dy[k] -= d_dHessenberg( k, i ) * d_dy[i];
        }

        d_dy[k] = d_dy[k] / d_dHessenberg( k, k );
    }
}

/****************************************************************
 *  Function to set the register the operator                    *
 ****************************************************************/
template<typename T>
void GMRESSolver<T>::registerOperator( std::shared_ptr<AMP::Operator::Operator> op )
{
    AMP_ASSERT( op );
    d_pOperator = op;
}

template<typename T>
void GMRESSolver<T>::resetOperator(
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

template<typename T>
void GMRESSolver<T>::computeInitialResidual( bool use_zero_guess,
                                             std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                                             std::shared_ptr<AMP::LinearAlgebra::Vector> u,
                                             std::shared_ptr<AMP::LinearAlgebra::Vector> tmp,
                                             std::shared_ptr<AMP::LinearAlgebra::Vector> res )
{
    if ( use_zero_guess ) {
        if ( d_bUsesPreconditioner && ( d_preconditioner_side == "left" ) ) {
            tmp->copyVector( f );
            d_pPreconditioner->apply( tmp, res );
        } else {
            res->copyVector( f );
        }
        u->setToScalar( static_cast<T>( 0.0 ) );
    } else {
        if ( d_bUsesPreconditioner && ( d_preconditioner_side == "left" ) ) {
            d_pOperator->residual( f, u, tmp );
            d_pPreconditioner->apply( tmp, res );
        } else {
            // this computes res = f - L(u)
            d_pOperator->residual( f, u, res );
        }
    }
}

template<typename T>
void GMRESSolver<T>::addCorrection( const int nr,
                                    std::shared_ptr<AMP::LinearAlgebra::Vector> z,
                                    std::shared_ptr<AMP::LinearAlgebra::Vector> z1,
                                    std::shared_ptr<AMP::LinearAlgebra::Vector> u )
{

    if ( d_bUsesPreconditioner && ( d_preconditioner_side == "right" ) ) {

        if ( !d_bFlexibleGMRES ) {
            z->setToScalar( static_cast<T>( 0.0 ) );

            for ( int i = 0; i <= nr; ++i ) {
                // this is V_m * y_m written as \sum_i i-th col(V_m) * i-th entry of y_m
                z->axpy( d_dy[i], *d_vBasis[i], *z );
            }

            // this solves M z1 = V_m * y_m (so z1= inv(M) * V_m * y_m)
            d_pPreconditioner->apply( z, z1 );
        } else {
            z1->setToScalar( static_cast<T>( 0.0 ) );
            for ( int i = 0; i <= nr; ++i ) {
                // this is Z_m * y_m written as \sum_i i-th col(Z_m) * i-th entry of y_m
                z1->axpy( d_dy[i], *d_zBasis[i], *z1 );
            }
        }
        u->axpy( 1.0, *z1, *u );

    } else {
        for ( int i = 0; i <= nr; ++i ) {
            u->axpy( d_dy[i], *d_vBasis[i], *u );
        }
    }
    u->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
}
} // namespace AMP::Solver
