#ifndef included_AMP_TrilinosNOXSolverParameters
#define included_AMP_TrilinosNOXSolverParameters

#include "AMP/solvers/SolverStrategyParameters.h"
DISABLE_WARNINGS
#include "AMP/solvers/trilinos/nox/PrePostOperator.h"
ENABLE_WARNINGS
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include <memory>


namespace AMP {
namespace Solver {

/**
 * Class TrilinosNOXSolverParameters provides a uniform mechanism to pass
 * initialization parameters to the PetscSNESSolver solver. It contains
 * shared pointers to a PertscKrylovSolver object and a vector
 * for initial guesses. All member variables are public.
 */
class TrilinosNOXSolverParameters : public SolverStrategyParameters
{
public:
    TrilinosNOXSolverParameters() {}
    explicit TrilinosNOXSolverParameters( std::shared_ptr<AMP::Database> db )
        : SolverStrategyParameters( db )
    {
    }
    virtual ~TrilinosNOXSolverParameters() {}

    AMP_MPI d_comm;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_pInitialGuess;
    std::shared_ptr<AMP::Operator::Operator> d_pLinearOperator;
    std::shared_ptr<AMP::Solver::SolverStrategy> d_preconditioner;
    std::shared_ptr<AMP::Solver::PrePostOperator> d_prePostOperator;

protected:
private:
};
} // namespace Solver
} // namespace AMP

#endif
