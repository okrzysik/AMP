#ifndef included_NonlinearKrylovAccelerator
#define included_NonlinearKrylovAccelerator


#include "AMP/utils/Database.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/Vector.h"
#include "NonlinearKrylovAcceleratorParameters.h"
#include "SolverStrategy.h"

namespace AMP {
namespace Solver {

/**
 * The NonlinearKrylovAccelerator class provides a C++ implementation of the NLKAIN (Nonlinear
 * Krylov Accelerator)
 * solver
 * as described and implemented by Carlson et. al.
 */
class NonlinearKrylovAccelerator : public SolverStrategy
{

public:
    /**
     * Main constructor.
     @param [in] params The parameters object of type NonlinearKrylovAcceleratorParameters contains
     a database object
     which must contain the
     following fields in addition to the fields expected by the base class SolverStrategy class:

       1. name:  max_vectors, type: integer, (required)
       acceptable values (non-negative integer values, typically 3 to 5)

       2. name:  angle_tolerance, type: double, (required)
       acceptable values (non-negative real values)

       3. name:  maximum_function_evals, type: integer, (required), default:
       acceptable values (non-negative integer values)

       4. name:  absolute_tolerance, type: double, (required), default: none
       acceptable values (non-negative real values)

       5. name:  relative_tolerance, type: double, (required), default: none
       acceptable values (non-negative real values)

       6. name:  freeze_pc, type: bool, (required), default: none
       acceptable values TRUE, FALSE)

     */
    explicit NonlinearKrylovAccelerator(
        AMP::shared_ptr<NonlinearKrylovAcceleratorParameters> params );

    /**
     * Default destructor
     */
    virtual ~NonlinearKrylovAccelerator();

    /**
     * Initialize the solution vector
     @param [in] parameters The parameters object (must be of type
     NonlinearKrylovAcceleratorParameters)
     Should contain a pointer to a solution guess vector.
     */
    void initialize( AMP::shared_ptr<SolverStrategyParameters> parameters ) override;

    /**
     * Provide the initial guess for the solver.
     * @param [in] initialGuess: shared pointer to the initial guess vector.
     */
    void setInitialGuess( AMP::shared_ptr<AMP::LinearAlgebra::Vector> initialGuess ) override;

    /**
     * This routine corrects the acceleration subspace using the vector f
     @param [in] f shared pointer to vector used to correct acceleration subspace
     */
    void correction( AMP::shared_ptr<AMP::LinearAlgebra::Vector> &f );

    /**
     * Resets the internal Krylov supsace, destroying all the internal vectors for a fresh start
     */
    void restart( void );

    /**
     * Should not be used at this time
     */
    void relax( void );

    /*!
     *  Get absolute tolerance for nonlinear solver.
     */
    double getAbsoluteTolerance() const;

    /*!
     *  Set absolute tolerance for nonlinear solver.
     */
    void setAbsoluteTolerance( double abs_tol );

    /*!
     *  Get relative tolerance for nonlinear solver.
     */
    double getRelativeTolerance() const;

    /*!
     *  Set relative tolerance for nonlinear solver.
     */
    void setRelativeTolerance( double rel_tol );
    /*!
     *  Get maximum iterations for nonlinear solver.
     */
    int getMaxNonlinearIterations() const;

    /*!
     *  Set maximum iterations for nonlinear solver.
     */
    void setMaxNonlinearIterations( int max_nli );

    /*!
     *  Get maximum function evaluations by nonlinear solver.
     */
    int getMaxFunctionEvaluations() const;

    /*!
     *  Set maximum function evaluations in nonlinear solver.
     */
    void setMaxFunctionEvaluations( int max_feval );

    /**
     * Solve the system \f$A(u) = f\f$.
     @param [in] f : shared pointer to right hand side vector
     @param [out] u : shared pointer to approximate computed solution
     */
    void solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                AMP::shared_ptr<AMP::LinearAlgebra::Vector> u ) override;

    /*!
     * Obtain number of nonlinear iterations.
     */
    int getNonlinearIterationCount() const { return d_iNonlinearIterationCount; }

    /**
     *
     */
    void putToDatabase( AMP::shared_ptr<AMP::Database> &db );

    /**
     * Sets the preconditioner that the solver will use.
     @param [in] preconditioner shared pointer to SolverStrategy object to use as preconditioner
     */
    void setPreconditioner( AMP::shared_ptr<AMP::Solver::SolverStrategy> preconditioner )
    {
        d_pPreconditioner = preconditioner;
    }

    /**
     * Return a shared pointer to the preconditioner being used.
     */
    AMP::shared_ptr<AMP::Solver::SolverStrategy> getPreconditioner( void )
    {
        return d_pPreconditioner;
    }

protected:
private:
    void getFromInput( const AMP::shared_ptr<AMP::Database> &db );

    bool d_bIsSubspace;                 /* boolean: a nonempty subspace */
    bool d_bContainsPendingVecs;        /* contains pending vectors -- boolean */
    int d_iMaximumNumberOfVectors;      /* maximum number of subspace vectors */
    double d_dVectorAngleDropTolerance; /* vector drop tolerance */

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_pvSolution;   /* correction vectors */
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_pvResidual;   /* correction vectors */
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_pvCorrection; /* correction vectors */

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> *d_pCorrectionVectors; /* correction vectors */
    AMP::shared_ptr<AMP::LinearAlgebra::Vector>
        *d_pFunctionDifferenceVectors;             /* function difference vectors */
    double **d_ppdFunctionDifferenceInnerProducts; /* matrix of w vector inner products */

    AMP::shared_ptr<AMP::Solver::SolverStrategy> d_pPreconditioner;

    /* Linked-list organization of the vector storage. */
    int d_iFirstVectorIndex; /* index of first subspace vector */
    int d_iLastVectorIndex;  /* index of last subspace vector */
    int d_iFreeVectorIndex;  /* index of the initial vector in free storage linked list */
    int *d_piNext;           /* next index link field */
    int *d_piPrevious;       /* previous index link field in doubly-linked subspace v */

    int d_iMaximumFunctionEvals;
    int d_iNonlinearIterationCount;

    double d_dAbsoluteTolerance;
    double d_dRelativeTolerance;

    bool d_bPrintResiduals;
    bool d_bSolverInitialized;
    bool d_bFreezePc;
};
} // namespace Solver
} // namespace AMP

#endif
