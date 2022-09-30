#ifndef included_AMP_GMRESSolver
#define included_AMP_GMRESSolver

#include "AMP/solvers/KrylovSolverParameters.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Array.h"

namespace AMP::Solver {

/**
 * The GMRESSolver class implements the GMRES method for non-symmetric linear systems
 * introduced by Saad and Schultz
 * Y. Saad and M.H. Schultz
 * "GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems",
 * SIAM J. Sci. Stat. Comput., 7:856-869, 1986.
 * doi:10.1137/0907058
 * If a preconditioner is provided right preconditioning is done
 */

template<typename T = double>
class GMRESSolver : public SolverStrategy
{
public:
    /**
     * default constructor
     */
    GMRESSolver();

    /**
     * main constructor
     @param [in] parameters The parameters object
     contains a database objects containing the following fields:

     1. type: double, name : relative_tolerance, default value of $1.0e-9$, relative tolerance for
    GMRES solver
    acceptable values (non-negative real values)

     2. type: bool, name : uses_preconditioner, default value false
        acceptable values (false, true),
        side effect: if false sets string pc_type to "none"

     3. type: string, name : pc_side, default value "RIGHT",
     acceptable values ("RIGHT", "LEFT" )
         active only when uses_preconditioner set to true
     */
    explicit GMRESSolver( std::shared_ptr<SolverStrategyParameters> parameters );

    /**
     * static create routine that is used by SolverFactory
     @param [in] parameters The parameters object
     contains a database objects with the fields listed for the constructor above
     */
    static std::unique_ptr<SolverStrategy>
    createSolver( std::shared_ptr<SolverStrategyParameters> solverStrategyParameters )
    {
        return std::make_unique<GMRESSolver<T>>( solverStrategyParameters );
    }

    /**
     * Default destructor
     */
    virtual ~GMRESSolver() = default;

    /**
     * Solve the system Au=f.
     * The implementation is based on the paper
     * "Implementations of the GMRES method", H. Walker
     * Computer Physics Communications 53 (1989) 311-320
     * @param [in] f : const shared pointer to right hand side vector
     * @param [out] u : shared pointer to approximate computed solution
     */
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                AMP::LinearAlgebra::Vector::shared_ptr u ) override;

    /**
     * Initialize the GMRESSolver. Should not be necessary for the user to call in general.
     * @param parameters
     */
    void initialize( std::shared_ptr<const SolverStrategyParameters> parameters ) override;

    /**
     * sets a shared pointer to a preconditioner object. The preconditioner is derived from
     * a SolverStrategy class
     * @param pc shared pointer to preconditioner
     */
    inline void setNestedSolver( std::shared_ptr<AMP::Solver::SolverStrategy> pc ) override
    {
        d_pPreconditioner = pc;
    }

    inline std::shared_ptr<AMP::Solver::SolverStrategy> getNestedSolver() override
    {
        return d_pPreconditioner;
    }

    /**
     * Register the operator that the solver will use during solves
     * @param [in] op shared pointer to operator $A()$ for equation \f$A(u) = f\f$
     */
    void registerOperator( std::shared_ptr<AMP::Operator::Operator> op ) override;

    /**
     * Resets the registered operator internally with new parameters if necessary
     * @param parameters    OperatorParameters object that is NULL by default
     */
    void
    resetOperator( std::shared_ptr<const AMP::Operator::OperatorParameters> parameters ) override;

protected:
    //! orthogonalize the vector against the existing vectors in the basis
    // stored internally. Store the coefficients of the Arnoldi
    // iteration internally in a upper Hessenberg matrix
    virtual void orthogonalize( const int k, std::shared_ptr<AMP::LinearAlgebra::Vector> v );

    //! apply the i-th Givens rotation to the k-th column of the Hessenberg matrix
    void applyGivensRotation( const int i, const int k );

    //! compute the Givens rotation required to zero out the sbub-diagonal element
    //! on the k-th column of the Hessenberg matrix and add it to the stored rotations
    //! Note that zero based indexing is being used
    void computeGivensRotation( const int k );

    //! perform a back solve for the upper triangular system generated for the least
    //! squares minimization problem
    //! @param[in] nr dimension of least squares system.
    void backwardSolve( const int nr );

    /**
     * update current approximation with the correction.
     * @param[in] nr dimension of least squares system.
     * @param[in] z shared pointer to temporary vector used if there is right preconditioning.
     * @param[in] z1 shared pointer to temporary vector used if there is right preconditioning.
     * @param[out] u shared pointer to approximate computed solution to correct.
     */
    void addCorrection( const int nr,
                        std::shared_ptr<AMP::LinearAlgebra::Vector> z,
                        std::shared_ptr<AMP::LinearAlgebra::Vector> z1,
                        std::shared_ptr<AMP::LinearAlgebra::Vector> u );

    /**
     * Compute initial residual for a GMRES cycle.
     * @param[in] use_zero_guess specify if zero should be used for the initial guess.
     * @param[in] f right hand side for residual.
     * @param[in,out] u initial guess if use_zero_guess is false (set to zero otherwise).
     * @param[in] tmp temporary vector needed if preconditiong is used.
     * @param[out] res residual vector.
     */
    void computeInitialResidual( bool use_zero_guess,
                                 std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                                 std::shared_ptr<AMP::LinearAlgebra::Vector> u,
                                 std::shared_ptr<AMP::LinearAlgebra::Vector> tmp,
                                 std::shared_ptr<AMP::LinearAlgebra::Vector> res );

    void getFromInput( std::shared_ptr<AMP::Database> db );

private:
    bool d_bRestart = false; //! whether to restart

    int d_iMaxKrylovDimension = 50; //! maximum dimension of the Krylov subspace before a restart or
                                    //! termination happens

    int d_restarts; //! logs number of times the solver is restarted

    //! string, determines orthogonalization method in Arnoldi
    //! options are "CGS", "MGS", "HR", where
    //! "CGS" : classical Gram-Schmidt ( fast but potentially unstable )
    //! "MGS" : modified Gram-Schmidt  ( stable )
    //! "HR" : Householder reflections (use when highly ill conditioned)
    std::string d_sOrthogonalizationMethod = "MGS";

    //! string, determining left, right or both side preconditioning
    //! this flag only applies if d_bUsesPreconditioner is true
    //! valid values are "left", "right", "both"
    //! currently only right and left are implemented
    std::string d_preconditioner_side = "right";

    //! boolean, for whether a preconditioner present or not
    bool d_bUsesPreconditioner = false;

    //! stores the upper Hessenberg matrix formed during the GMRES iteration
    AMP::Array<T> d_dHessenberg;

    //! stores the cosine and sine values required for Givens rotations
    // when Givens is being used
    std::vector<T> d_dcos;
    std::vector<T> d_dsin;

    //! stores the right hand side for the Hessenberg least squares system
    std::vector<T> d_dw;

    //! stores the solution for the least squares system
    std::vector<T> d_dy;

    //! shared pointer to preconditioner if it exists
    std::shared_ptr<AMP::Solver::SolverStrategy> d_pPreconditioner;

    //! stores the orthonormal basis for the Krylov space
    //! we do not preallocate by default
    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_vBasis;
};
} // namespace AMP::Solver

#endif
