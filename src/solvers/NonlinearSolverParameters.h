#ifndef included_NonlinearSolverParameters_h_
#define included_NonlinearSolverParameters_h_

#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/utils/AMP_MPI.h"

namespace AMP::LinearAlgebra {
class Vector;
}

namespace AMP::Solver {

class SolverStrategy;

class NonlinearSolverParameters : public SolverStrategyParameters
{
public:
    NonlinearSolverParameters() = default;

    /**
     * Construct and initialize a parameter list according to input
     * data.  Guess what the required and optional keywords are.
     */
    explicit NonlinearSolverParameters( std::shared_ptr<AMP::Database> db )
        : SolverStrategyParameters( db )
    {
    }

    ~NonlinearSolverParameters() = default;

    AMP_MPI d_comm;

    //! pointer to linear solver, default of nullptr
    //    std::shared_ptr<AMP::Solver::SolverStrategy> d_linear_solver;
    //! pointer to preconditioner, default of nullptr
    std::shared_ptr<AMP::Solver::SolverStrategy> d_pNestedSolver;

    //! initial guess for solver
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_pInitialGuess;
};
} // namespace AMP::Solver

#endif
