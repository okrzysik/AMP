#ifndef included_AMP_DiagonalSolver
#define included_AMP_DiagonalSolver

#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/utils/AMP_MPI.h"

#include <string>

namespace AMP::Solver {

/**
 * The DiagonalSolver class simply does diagonal scaling
 */

template<typename T = double>
class DiagonalSolver : public SolverStrategy
{
public:
    /**
     * default constructor
     */
    DiagonalSolver() = default;

    /**
     * main constructor
     @param [in] params The parameters object
     contains a database objects containing the following fields:

     1. type: double, name : relative_tolerance, default value of $1.0e-9$, relative tolerance for
    CG solver
    acceptable values (non-negative real values)

     2. type: bool, name : uses_nested_solver, default value false
        acceptable values (false, true),
        side effect: if false sets string nested_solver_type to "none"
     */
    explicit DiagonalSolver( std::shared_ptr<SolverStrategyParameters> params );

    /**
     * static create routine that is used by SolverFactory
     @param [in] params The parameters object
     contains a database objects with the fields listed for the constructor above
     */
    static std::unique_ptr<SolverStrategy>
    createSolver( std::shared_ptr<SolverStrategyParameters> params )
    {
        return std::make_unique<DiagonalSolver<T>>( params );
    }

    /**
     * Default destructor
     */
    virtual ~DiagonalSolver() = default;

    std::string type() const override { return "DiagonalSolver"; }

    /**
     * Solve the system \f$Au = 0\f$.
     * @param [in] f : shared pointer to right hand side vector
     * @param [out] u : shared pointer to approximate computed solution
     */
    void apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                std::shared_ptr<AMP::LinearAlgebra::Vector> u ) override;

    /**
     * Initialize the DiagonalSolver. Should not be necessary for the user to call in general.
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
     * @param [in] op shared pointer to operator \f$A()\f$ for equation \f$A(u) = f\f$
     */
    void registerOperator( std::shared_ptr<AMP::Operator::Operator> op ) override;

    /**
     * Resets the registered operator internally with new parameters if necessary
     * @param params    OperatorParameters object that is NULL by default
     */
    void resetOperator( std::shared_ptr<const AMP::Operator::OperatorParameters> params ) override;

protected:
    void getFromInput( std::shared_ptr<AMP::Database> db );

private:
    bool d_bUsesNestedSolver = false;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_pDiagonalInverse;
    std::shared_ptr<AMP::Solver::SolverStrategy> d_pNestedSolver;
};
} // namespace AMP::Solver

#endif
