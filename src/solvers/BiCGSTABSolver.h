#ifndef included_AMP_BiCGSTABSolver
#define included_AMP_BiCGSTABSolver

#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/utils/AMP_MPI.h"


namespace AMP::Solver {


/**
 * The BiCGSTABSolver class implements the BiCGSTAB method for non-symmetric linear systems
 * introduced by H. van der Vorst
 *    Van der Vorst, H. A., "Bi-CGSTAB: A Fast and Smoothly Converging Variant of
 *    Bi-CG for the Solution of Nonsymmetric Linear Systems". SIAM J. Sci. and Stat. Comput.
 *    13 (2): 631-644 (1992). doi:10.1137/0913035.
 * If a preconditioner is provided right preconditioning is done
 */
template<typename T = double>
class BiCGSTABSolver : public SolverStrategy
{
public:
    /**
     * default constructor
     */
    BiCGSTABSolver();

    /**
     * main constructor
     @param [in] params The parameters object
     contains a database objects containing the following fields:

     1. type: double, name : relative_tolerance, default value of $1.0e-9$, relative tolerance for
    CG solver
    acceptable values (non-negative real values)

     2. type: bool, name : uses_preconditioner, default value false
        acceptable values (false, true),
        side effect: if false sets string pc_type to "none"

     3. type: string, name : pc_side, default value "RIGHT",
     acceptable values ("RIGHT", "LEFT", "SYMMETRIC" )
         active only when uses_preconditioner set to true
     */
    explicit BiCGSTABSolver( std::shared_ptr<SolverStrategyParameters> params );

    /**
     * static create routine that is used by SolverFactory
     @param [in] params The parameters object
     contains a database objects with the fields listed for the constructor above
     */
    static std::unique_ptr<SolverStrategy>
    createSolver( std::shared_ptr<SolverStrategyParameters> params )
    {
        return std::make_unique<BiCGSTABSolver<T>>( params );
    }

    /**
     * Default destructor
     */
    virtual ~BiCGSTABSolver();

    std::string type() const override { return "BiCGSTABSolver"; }

    /**
     * Solve the system \f$Au = 0\f$.
     * @param [in] f : shared pointer to right hand side vector
     * @param [out] u : shared pointer to approximate computed solution
     */
    void apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                std::shared_ptr<AMP::LinearAlgebra::Vector> u ) override;

    /**
     * Initialize the BiCGSTABSolver. Should not be necessary for the user to call in general.
     * @param params
     */
    void initialize( std::shared_ptr<const SolverStrategyParameters> params ) override;

    /**
     * sets a shared pointer to a preconditioner object. The preconditioner is derived from
     * a SolverStrategy class
     * @param pc shared pointer to preconditioner
     */
    inline void setNestedSolver( std::shared_ptr<AMP::Solver::SolverStrategy> pc ) override
    {
        d_pNestedSolver = pc;
    }

    inline std::shared_ptr<AMP::Solver::SolverStrategy> getNestedSolver() override
    {
        return d_pNestedSolver;
    }

    /**
     * Register the operator that the solver will use during solves
     * @param [in] op shared pointer to operator $A()$ for equation \f$A(u) = f\f$
     */
    void registerOperator( std::shared_ptr<AMP::Operator::Operator> op ) override;

    /**
     * Resets the registered operator internally with new parameters if necessary
     * @param params    OperatorParameters object that is NULL by default
     */
    void resetOperator( std::shared_ptr<const AMP::Operator::OperatorParameters> params ) override;

protected:
    void getFromInput( std::shared_ptr<AMP::Database> db );
    void allocateScratchVectors( std::shared_ptr<const AMP::LinearAlgebra::Vector> u );

private:
    int d_restarts = 0; //! number of times the solver is restarted

    bool d_bUsesPreconditioner = false;

    //! scratch vectors required for BiCGSTAB
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_r;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_r_tilde;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_p;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_p_hat;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_s;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_s_hat;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_t;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_v;

    std::shared_ptr<AMP::Solver::SolverStrategy> d_pNestedSolver;
};
} // namespace AMP::Solver

#endif
