#ifndef included_AMP_CGSolver
#define included_AMP_CGSolver

#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/utils/AMP_MPI.h"

#include <string>

namespace AMP::Solver {

/**
 * The CGSolver class implements the Conjugate Gradient method namely the 2-term recurrence variant.
 * M.R. Hestenes, E. Stiefel. "Methods of conjugate gradients for solving linear systems"
 * J. Res. Natl. Bur. Stand., 49 (1952), pp. 409-436
 *
 * In addition it implements the IPCG variant  developed in
 * Golub, Gene H.; Ye, Qiang (1999). "Inexact Preconditioned Conjugate Gradient Method with
 * Inner-Outer Iteration". SIAM Journal on Scientific Computing 21 (4): 1305.
 * doi:10.1137/S1064827597323415 (http://dx.doi.org/10.1137%2FS1064827597323415) .
 */

template<typename T = double>
class CGSolver : public SolverStrategy
{
public:
    /**
     * default constructor
     */
    CGSolver() = default;

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
     */
    explicit CGSolver( std::shared_ptr<SolverStrategyParameters> params );

    /**
     * static create routine that is used by SolverFactory
     @param [in] params The parameters object
     contains a database objects with the fields listed for the constructor above
     */
    static std::unique_ptr<SolverStrategy>
    createSolver( std::shared_ptr<SolverStrategyParameters> params )
    {
        return std::make_unique<CGSolver<T>>( params );
    }

    /**
     * Default destructor
     */
    virtual ~CGSolver() = default;

    std::string type() const override { return "CGSolver"; }

    /**
     * Solve the system \f$Au = 0\f$.
     * @param [in] f : shared pointer to right hand side vector
     * @param [out] u : shared pointer to approximate computed solution
     */
    void apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                std::shared_ptr<AMP::LinearAlgebra::Vector> u ) override;

    /**
     * Initialize the CGSolver. Should not be necessary for the user to call in general.
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
        d_pPreconditioner = pc;
    }

    inline std::shared_ptr<AMP::Solver::SolverStrategy> getNestedSolver() override
    {
        return d_pPreconditioner;
    }

    /**
     * Resets the registered operator internally with new parameters if necessary
     * @param params    OperatorParameters object that is NULL by default
     */
    void resetOperator( std::shared_ptr<const AMP::Operator::OperatorParameters> params ) override;

protected:
    void getFromInput( std::shared_ptr<AMP::Database> db );

private:
    T d_dDivergenceTolerance = 1e3;

    bool d_bUsesPreconditioner = false;

    //! use flexible CG if true
    bool d_bFlexibleCG = false;

    //! maximum dimension of the stored search space for FCG
    int d_max_dimension = 0;

    //! variant being used, can be one of "pcg", "ipcg", or "fcg"
    std::string d_sVariant = "pcg";

    std::shared_ptr<AMP::Solver::SolverStrategy> d_pPreconditioner;

    //! stores the search directions for IPCG/FCG if needed
    //! we do not preallocate by default
    std::vector<std::shared_ptr<AMP::LinearAlgebra::Vector>> d_vDirs;
};
} // namespace AMP::Solver

#endif
