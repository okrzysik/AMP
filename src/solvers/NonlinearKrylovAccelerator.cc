#include "NonlinearKrylovAccelerator.h"

#include "utils/Utilities.h"

#include "utils/Utilities.h"

#define EOL -1

namespace AMP {
namespace Solver {

NonlinearKrylovAccelerator::NonlinearKrylovAccelerator(
    AMP::shared_ptr<NonlinearKrylovAcceleratorParameters> params )
    : SolverStrategy( params )
{
    int j, n;

    d_pCorrectionVectors         = nullptr;
    d_pFunctionDifferenceVectors = nullptr;
    d_iNonlinearIterationCount   = 0;
    d_iMaximumFunctionEvals      = 50;
    d_dAbsoluteTolerance         = 1.0e-12;
    d_dRelativeTolerance         = 1.0e-8;

    d_bPrintResiduals    = false;
    d_bSolverInitialized = false;
    d_bFreezePc          = true;

    getFromInput( params->d_db );

    d_pPreconditioner = params->d_pPreconditioner;

    n = d_iMaximumNumberOfVectors + 1;

    d_ppdFunctionDifferenceInnerProducts    = new double *[n];
    d_ppdFunctionDifferenceInnerProducts[0] = new double[n * n];

    for ( j = 1; j < n; j++ ) {
        d_ppdFunctionDifferenceInnerProducts[j] = d_ppdFunctionDifferenceInnerProducts[j - 1] + n;
    }

    d_piNext     = new int[n];
    d_piPrevious = new int[n];

    if ( params->d_pInitialGuess.get() != nullptr ) {
        initialize( params );
    }

    restart();
}


NonlinearKrylovAccelerator::~NonlinearKrylovAccelerator( void )
{
    if ( d_ppdFunctionDifferenceInnerProducts != nullptr ) {
        delete d_ppdFunctionDifferenceInnerProducts[0];
        delete[] d_ppdFunctionDifferenceInnerProducts;
    }

    delete[] d_piNext;
    delete[] d_piPrevious;

    /* freeVectorComponents is no longer necessary
    for (int j = 0; j < d_iMaximumNumberOfVectors+1; j++)
    {
        if( (d_pCorrectionVectors!=NULL)&&(d_pCorrectionVectors[j].get()!=NULL) )
        {
            d_pCorrectionVectors[j]->freeVectorComponents();
        }

        if((d_pFunctionDifferenceVectors!=NULL) && (d_pFunctionDifferenceVectors[j].get()!=NULL))
        {
            d_pFunctionDifferenceVectors[j]->freeVectorComponents();
        }
    }*/

    delete[] d_pCorrectionVectors;
    delete[] d_pFunctionDifferenceVectors;
}


void NonlinearKrylovAccelerator::getFromInput( const AMP::shared_ptr<AMP::Database> &db )
{
    if ( db->keyExists( "max_vectors" ) ) {
        d_iMaximumNumberOfVectors = db->getInteger( "max_vectors" );
    } else {
        AMP_ERROR( "NonlinearKrylovAccelerator"
                   << " -- Key data `max_vectors'"
                   << " missing in input." );
    }

    if ( db->keyExists( "angle_tolerance" ) ) {
        d_dVectorAngleDropTolerance = db->getDouble( "angle_tolerance" );
    } else {
        AMP_ERROR( "NonlinearKrylovAccelerator"
                   << " -- Key data `angle_tolerance'"
                   << " missing in input." );
    }

    if ( db->keyExists( "maximum_function_evals" ) ) {
        d_iMaximumFunctionEvals = db->getInteger( "maximum_function_evals" );
    }

    if ( db->keyExists( "absolute_tolerance" ) ) {
        d_dAbsoluteTolerance = db->getDouble( "absolute_tolerance" );
    }

    if ( db->keyExists( "relative_tolerance" ) ) {
        d_dRelativeTolerance = db->getDouble( "relative_tolerance" );
    }

    if ( db->keyExists( "freeze_pc" ) ) {
        d_bFreezePc = db->getBool( "freeze_pc" );
    }

    if ( d_iDebugPrintInfoLevel > 0 ) {
        d_bPrintResiduals = true;
    }

    AMP_INSIST( d_iMaximumNumberOfVectors > 0, "The maximum number of vectors must be positive" );
    AMP_INSIST( d_dVectorAngleDropTolerance > 0.0, "The tolerance in angle must be positive" );
}


void NonlinearKrylovAccelerator::setInitialGuess(
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> initialGuess )
{
    size_t n;

    d_pvSolution = initialGuess;

    n = d_iMaximumNumberOfVectors + 1;

    if ( !d_bSolverInitialized ) {
        d_pCorrectionVectors = new AMP::shared_ptr<AMP::LinearAlgebra::Vector>[n];

        for ( size_t j = 0; j < n; j++ ) {
            d_pCorrectionVectors[j] = d_pvSolution->cloneVector();
        }

        d_pFunctionDifferenceVectors = new AMP::shared_ptr<AMP::LinearAlgebra::Vector>[n];

        for ( size_t j = 0; j < n; j++ ) {
            d_pFunctionDifferenceVectors[j] = d_pvSolution->cloneVector();
        }

        d_pvResidual   = d_pvSolution->cloneVector();
        d_pvCorrection = d_pvSolution->cloneVector();

        d_bSolverInitialized = true;
    }

    d_pvResidual->setToScalar( 0.0 );
    d_pvCorrection->setToScalar( 0.0 );
}


void NonlinearKrylovAccelerator::initialize( AMP::shared_ptr<SolverStrategyParameters> parameters )
{
    AMP::shared_ptr<NonlinearKrylovAcceleratorParameters> params =
        AMP::dynamic_pointer_cast<NonlinearKrylovAcceleratorParameters>( parameters );
    setInitialGuess( params->d_pInitialGuess );
}


void NonlinearKrylovAccelerator::correction( AMP::shared_ptr<AMP::LinearAlgebra::Vector> &f )
{
    int i, j, k, new_loc;
    double s;
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> v, w;
    double *hj, *c;
    /*
     *  UPDATE THE ACCELERATION SUBSPACE
     */

    if ( d_bContainsPendingVecs ) {

        /* next function difference w_1 */
        w = d_pFunctionDifferenceVectors[d_iFirstVectorIndex];
        w->axpy( -1.0, *f, *w );
        s = w->L2Norm();

        /* If the function difference is 0, we can't update the subspace with
        this data; so we toss it out and continue.  In this case it is likely
        that the outer iterative solution procedure has gone badly awry
        (unless the function value is itself 0), and we merely want to do
        something reasonable here and hope that situation is detected on the
        outside. */
        if ( s == 0.0 ) {
            relax();
        }
    }

    if ( d_bContainsPendingVecs ) {

        /* Normalize w_1 and apply same factor to v_1. */
        v = d_pCorrectionVectors[d_iFirstVectorIndex];

        v->scale( 1.0 / s, *v );
        w->scale( 1.0 / s, *w );

        /* Update H. */
        for ( k = d_piNext[d_iFirstVectorIndex]; k != EOL; k = d_piNext[k] ) {
            d_ppdFunctionDifferenceInnerProducts[d_iFirstVectorIndex][k] =
                w->dot( *d_pFunctionDifferenceVectors[k] );
        }

        /*
         *  CHOLESKI FACTORIZATION OF H = W^t W
         *  original matrix kept in the upper triangle (implicit unit diagonal)
         *  lower triangle holds the factorization
        */

        /* Trivial initial factorization stage. */
        int nvec                                                                       = 1;
        d_ppdFunctionDifferenceInnerProducts[d_iFirstVectorIndex][d_iFirstVectorIndex] = 1.0;

        for ( k = d_piNext[d_iFirstVectorIndex]; k != EOL; k = d_piNext[k] ) {

            /* Maintain at most MVEC vectors. */
            if ( ++nvec > d_iMaximumNumberOfVectors ) {
                /* Drop the last vector and update the free storage list. */
                AMP_INSIST( d_iLastVectorIndex == k, "d_iLastVectorIndex not equal to k" );
                d_piNext[d_iLastVectorIndex] = d_iFreeVectorIndex;
                d_iFreeVectorIndex           = k;
                d_iLastVectorIndex           = d_piPrevious[k];
                d_piNext[d_iLastVectorIndex] = EOL;
                break;
            }

            /* Single stage of Choleski factorization. */
            double *hk = d_ppdFunctionDifferenceInnerProducts[k]; /* row k of H */
            double hkk = 1.0;
            for ( j = d_iFirstVectorIndex; j != k; j = d_piNext[j] ) {
                hj         = d_ppdFunctionDifferenceInnerProducts[j]; /* row j of H */
                double hkj = hj[k];
                for ( i = d_iFirstVectorIndex; i != j; i = d_piNext[i] )
                    hkj -= hk[i] * hj[i];
                hkj /= hj[j];
                hk[j] = hkj;
                hkk -= hkj * hkj;
            }

            if ( hkk > pow( d_dVectorAngleDropTolerance, 2 ) ) {
                hk[k] = sqrt( hkk );
            } else {
                /* The current w nearly lies in the span of the previous vectors: */
                /* Drop this vector, */
                AMP_INSIST( d_piPrevious[k] != EOL, "The previous vector index equal to EOL" );
                d_piNext[d_piPrevious[k]] = d_piNext[k];
                if ( d_piNext[k] == EOL )
                    d_iLastVectorIndex = d_piPrevious[k];
                else
                    d_piPrevious[d_piNext[k]] = d_piPrevious[k];
                /* update the free storage list, */
                d_piNext[k]        = d_iFreeVectorIndex;
                d_iFreeVectorIndex = k;
                /* back-up and move on to the piNext vector. */
                k = d_piPrevious[k];
                nvec--;
            }
        }

        AMP_INSIST( d_iFirstVectorIndex != EOL, "IFirstVectorIndex vector index equal to EOL" );
        d_bIsSubspace = true; /* the acceleration subspace isn't empty */
    }

    /*
     *  ACCELERATED CORRECTION
     */

    /* Locate storage for the new vectors. */
    AMP_INSIST( d_iFreeVectorIndex != EOL, "d_iFreeVectorIndex is equal to EOL" );
    new_loc            = d_iFreeVectorIndex;
    d_iFreeVectorIndex = d_piNext[d_iFreeVectorIndex];

    /* Save the original f for the next call. */
    d_pFunctionDifferenceVectors[new_loc]->copyVector( f );

    if ( d_bIsSubspace ) {
        c = new double[( d_iMaximumNumberOfVectors + 1 )];

        AMP_INSIST( c != nullptr, "c is NULL" );

        /* Project f onto the span of the w vectors: */
        /* forward substitution */
        for ( j = d_iFirstVectorIndex; j != EOL; j = d_piNext[j] ) {
            double cj = f->dot( *d_pFunctionDifferenceVectors[j] );

            for ( i = d_iFirstVectorIndex; i != j; i = d_piNext[i] ) {
                cj -= d_ppdFunctionDifferenceInnerProducts[j][i] * c[i];
            }

            c[j] = cj / d_ppdFunctionDifferenceInnerProducts[j][j];
        }

        /* backward substitution */
        for ( j = d_iLastVectorIndex; j != EOL; j = d_piPrevious[j] ) {
            double cj = c[j];
            for ( i = d_iLastVectorIndex; i != j; i = d_piPrevious[i] )
                cj -= d_ppdFunctionDifferenceInnerProducts[i][j] * c[i];
            c[j] = cj / d_ppdFunctionDifferenceInnerProducts[j][j];
        }

        /* The accelerated correction */
        for ( k = d_iFirstVectorIndex; k != EOL; k = d_piNext[k] ) {
            f->axpy( c[k], *d_pCorrectionVectors[k], *f );
            f->axpy( -c[k], *d_pFunctionDifferenceVectors[k], *f );
        }

        delete[] c;
    }

    /* Save the accelerated correction for the next call. */
    d_pCorrectionVectors[new_loc]->copyVector( f );

    /* Prepend the new vectors to the list. */
    d_piPrevious[new_loc] = EOL;
    d_piNext[new_loc]     = d_iFirstVectorIndex;

    if ( d_iFirstVectorIndex == EOL ) {
        d_iLastVectorIndex = new_loc;
    } else {
        d_piPrevious[d_iFirstVectorIndex] = new_loc;
    }

    d_iFirstVectorIndex = new_loc;

    /* The original f and accelerated correction are cached for the next call. */
    d_bContainsPendingVecs = true;
}


void NonlinearKrylovAccelerator::solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                                        AMP::shared_ptr<AMP::LinearAlgebra::Vector>
                                            u )
{
    AMP_INSIST( d_pOperator != nullptr, "Operator cannot be NULL" );
    AMP_INSIST( d_pPreconditioner.get() != nullptr, "Preconditioning operator cannot be NULL" );

    double residual_norm = 1.0e10;

    d_iNonlinearIterationCount = 0;

    // make the internal solution vector point to the output solution so we don't have to make a
    // copy
    d_pvSolution = u;

    // compute residual
    d_pOperator->residual( f, d_pvSolution, d_pvResidual );
    residual_norm = d_pvResidual->L2Norm();

    if ( d_bPrintResiduals ) {
        AMP::pout << "NonlinearKrylovAccelerator::solve: iteration : " << d_iNonlinearIterationCount
                  << ", residual: " << residual_norm << std::endl;
    }

    AMP::shared_ptr<AMP::Operator::OperatorParameters> pc_parameters =
        d_pOperator->getParameters( "Jacobian", d_pvSolution );
    AMP::shared_ptr<AMP::Operator::Operator> pc_operator = d_pPreconditioner->getOperator();

    AMP_INSIST( pc_operator.get() != nullptr,
                "NonlinearKrylovAccelerator::solve: preconditioning operator cannot be NULL" );

    // if using a frozen preconditioner set it up iFirstVectorIndex
    if ( d_bFreezePc ) {
        d_pPreconditioner->resetOperator( pc_parameters );
    }

    while ( ( d_iNonlinearIterationCount < d_iMaxIterations ) &&
            ( residual_norm > d_dAbsoluteTolerance ) ) {
        if ( !d_bFreezePc ) {
            pc_parameters = d_pOperator->getParameters( "Jacobian", d_pvSolution );
            d_pPreconditioner->resetOperator( pc_parameters );
        }

        // apply the preconditioner
        d_pPreconditioner->solve( d_pvResidual, d_pvCorrection );

        if ( d_iDebugPrintInfoLevel > 3 ) {
            // compute residual
            d_pOperator->residual( f, d_pvSolution, d_pvResidual );
            residual_norm = d_pvResidual->L2Norm();
            std::cout << "NonlinearKrylovAccelerator::solve: L2 norm of preconditioner residual "
                      << residual_norm << std::endl;
        }


        // compute AIN correction
        this->correction( d_pvCorrection );

        if ( d_iDebugPrintInfoLevel > 3 ) {
            double pcSolutionNorm = d_pvCorrection->L2Norm();
            std::cout << "NonlinearKrylovAccelerator::solve: L2 norm of correction from NLKAIN "
                      << pcSolutionNorm << std::endl;
        }

        // correct current solution
        d_pvSolution->axpy( 1.0, *d_pvCorrection, *d_pvSolution );

        // compute the residual
        d_pOperator->residual( f, d_pvSolution, d_pvResidual );

        if ( d_iDebugPrintInfoLevel > 3 ) {
            double pcSolutionNorm = d_pvResidual->L2Norm();
            std::cout << "NonlinearKrylovAccelerator::solve: L2 norm of corrected residual "
                      << pcSolutionNorm << std::endl;
        }

        residual_norm = d_pvResidual->L2Norm();
        d_iNonlinearIterationCount++;

        if ( d_bPrintResiduals ) {
            AMP::pout << "Nonlinear Krylov iteration : " << d_iNonlinearIterationCount
                      << ", residual: " << residual_norm << std::endl;
        }
    }
}


void NonlinearKrylovAccelerator::restart( void )
{
    int k;

    /* No vectors are stored. */
    d_iFirstVectorIndex    = EOL;
    d_iLastVectorIndex     = EOL;
    d_bIsSubspace          = false;
    d_bContainsPendingVecs = false;

    /* Initialize the free storage linked list. */
    d_iFreeVectorIndex = 0;

    for ( k = 0; k < d_iMaximumNumberOfVectors; k++ ) {
        d_piNext[k] = k + 1;
    }

    d_piNext[d_iMaximumNumberOfVectors] = EOL;

    d_iNonlinearIterationCount = 0;
}


void NonlinearKrylovAccelerator::relax( void )
{
    if ( d_bContainsPendingVecs ) {
        /* Drop the initial slot where the pending vectors are stored. */
        AMP_INSIST( d_iFirstVectorIndex >= 0, "d_iFirstVectorIndex is not positive" );
        int new_loc         = d_iFirstVectorIndex;
        d_iFirstVectorIndex = d_piNext[d_iFirstVectorIndex];
        if ( d_iFirstVectorIndex == EOL ) {
            d_iLastVectorIndex = EOL;
        } else {
            d_piPrevious[d_iFirstVectorIndex] = EOL;
        }

        /* Update the free storage list. */
        d_piNext[new_loc]      = d_iFreeVectorIndex;
        d_iFreeVectorIndex     = new_loc;
        d_bContainsPendingVecs = false;
    }
}


/*
*************************************************************************
*                                                                       *
*  Access functions for parameters that control solver behavior.        *
*                                                                       *
*************************************************************************
*/
double NonlinearKrylovAccelerator::getAbsoluteTolerance() const { return ( d_dAbsoluteTolerance ); }

void NonlinearKrylovAccelerator::setAbsoluteTolerance( double abs_tol )
{
    d_dAbsoluteTolerance = abs_tol;
}

double NonlinearKrylovAccelerator::getRelativeTolerance() const { return ( d_dRelativeTolerance ); }

void NonlinearKrylovAccelerator::setRelativeTolerance( double rel_tol )
{
    d_dRelativeTolerance = rel_tol;
}


int NonlinearKrylovAccelerator::getMaxNonlinearIterations() const { return ( d_iMaxIterations ); }

void NonlinearKrylovAccelerator::setMaxNonlinearIterations( int max_nli )
{
    d_iMaxIterations = max_nli;
}

int NonlinearKrylovAccelerator::getMaxFunctionEvaluations() const
{
    return ( d_iMaximumFunctionEvals );
}

void NonlinearKrylovAccelerator::setMaxFunctionEvaluations( int max_feval )
{
    d_iMaximumFunctionEvals = max_feval;
}


void NonlinearKrylovAccelerator::putToDatabase( AMP::shared_ptr<AMP::Database> &db )
{

    AMP_INSIST( db.get() != nullptr, "database object cannot be NULL" );

    db->putInteger( "d_MaxIterations", d_iMaxIterations );
    db->putInteger( "d_iMaximumFunctionEvals", d_iMaximumFunctionEvals );

    db->putDouble( "d_dAbsoluteTolerance", d_dAbsoluteTolerance );
    db->putDouble( "d_dRelativeTolerance", d_dRelativeTolerance );
}
}
}
