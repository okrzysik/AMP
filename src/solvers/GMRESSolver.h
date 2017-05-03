#ifndef included_AMP_GMRESSolver
#define included_AMP_GMRESSolver

#include "solvers/KrylovSolverParameters.h"
#include "solvers/SolverStrategy.h"
#include "utils/AMP_MPI.h"
#include "utils/Array.h"

namespace AMP {
namespace Solver {

/**
 * The GMRESSolver class implements the GMRES method for non-symmetric linear systems
 * introduced by Saad and Schultz
 * Y. Saad and M.H. Schultz
 * "GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems",
 * SIAM J. Sci. Stat. Comput., 7:856-869, 1986.
 * doi:10.1137/0907058
 * If a preconditioner is provided right preconditioning is done
 */

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
    explicit GMRESSolver( AMP::shared_ptr<KrylovSolverParameters> parameters );

    /**
     * Default destructor
     */
    virtual ~GMRESSolver();

    /**
     * Solve the system Au=f.
     * The implementation is based on the paper
     * "Implementations of the GMRES method", H. Walker
     * Computer Physics Communications 53 (1989) 311-320
     * @param [in] f : const shared pointer to right hand side vector
     * @param [out] u : shared pointer to approximate computed solution
     */
    void solve( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                AMP::LinearAlgebra::Vector::shared_ptr u ) override;

    /**
     * Initialize the GMRESSolver. Should not be necessary for the user to call in general.
     * @param parameters
     */
    void initialize( AMP::shared_ptr<SolverStrategyParameters> const parameters ) override;

    /**
     * returns a shared pointer to a preconditioner object. The preconditioner is derived from
     * a SolverStrategy class
     */
    inline AMP::shared_ptr<AMP::Solver::SolverStrategy> getPreconditioner( void )
    {
        return d_pPreconditioner;
    }

    /**
     * sets a shared pointer to a preconditioner object. The preconditioner is derived from
     * a SolverStrategy class
     * @param pc shared pointer to preconditioner
     */
    inline void setPreconditioner( AMP::shared_ptr<AMP::Solver::SolverStrategy> pc )
    {
        d_pPreconditioner = pc;
    }

    /**
     * Register the operator that the solver will use during solves
     * @param [in] op shared pointer to operator $A()$ for equation \f$A(u) = f\f$
     */
    void registerOperator( const AMP::shared_ptr<AMP::Operator::Operator> op ) override;

    /**
     * Resets the registered operator internally with new parameters if necessary
     * @param parameters    OperatorParameters object that is NULL by default
     */
    void resetOperator( const AMP::shared_ptr<AMP::Operator::OperatorParameters> parameters ) override;

protected:
    //! orthogonalize the vector against the existing vectors in the basis
    // stored internally. Store the coefficients of the Arnoldi
    // iteration internally in a upper Hessenberg matrix
    virtual void orthogonalize( AMP::shared_ptr<AMP::LinearAlgebra::Vector> v );

    //! apply the i-th Givens rotation to the k-th column of the Hessenberg matrix
    void applyGivensRotation( const int i, const int k );

    //! compute the Givens rotation required to zero out the sbub-diagonal element
    //! on the k-th column of the Hessenberg matrix and add it to the stored rotations
    //! Note that zero based indexing is being used
    void computeGivensRotation( const int k );

    //! perform a back solve for the upper triangular system generated for the least
    //! squares minimization problem
    void backwardSolve( void );

    void getFromInput( const AMP::shared_ptr<AMP::Database> &db );

private:
    AMP_MPI d_comm;

    double d_dRelativeTolerance; //! relative tolerance to converge to

    bool d_bRestart; //! whether to restart

    int d_iMaxKrylovDimension; //! maximum dimension of the Krylov subspace before a restart or
                               //! termination happens

    int d_restarts; //! logs number of times the solver is restarted

    int d_nr; // dimension of the least squares system

    //! string, determines orthogonalization method in Arnoldi
    //! options are "CGS", "MGS", "HR", where
    //! "CGS" : classical Gram-Schmidt ( fast but potentially unstable )
    //! "MGS" : modified Gram-Schmidt  ( stable )
    //! "HR" : Householder reflections (use when highly ill conditioned)
    std::string d_sOrthogonalizationMethod;

    //! boolean, for whether a preconditioner present or not
    bool d_bUsesPreconditioner;

    //! stores the upper Hessenberg matrix formed during the GMRES iteration
    AMP::Array<double> d_dHessenberg;

    //! stores the cosine and sine values required for Givens rotations
    // when Givens is being used
    std::vector<double> d_dcos;
    std::vector<double> d_dsin;

    //! stores the right hand side for the Hessenberg least squares system
    std::vector<double> d_dw;

    //! stores the solution for the least squares system
    std::vector<double> d_dy;

    //! shared pointer to preconditioner if it exists
    AMP::shared_ptr<AMP::Solver::SolverStrategy> d_pPreconditioner;

    //! stores the orthonormal basis for the Krylov space
    //! we do not preallocate by default
    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_vBasis;
};
}
}

#endif
