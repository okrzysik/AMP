#ifndef included_AMP_TrilinosNOXSolverParameters
#define included_AMP_TrilinosNOXSolverParameters

#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/solvers/trilinos/nox/PrePostOperator.h"
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
    AMP::LinearAlgebra::Vector::shared_ptr d_pInitialGuess;
    AMP::Operator::Operator::shared_ptr d_pLinearOperator;
    AMP::Solver::SolverStrategy::shared_ptr d_preconditioner;
    AMP::Solver::PrePostOperator::shared_ptr d_prePostOperator;

protected:
private:
};
} // namespace Solver
} // namespace AMP

#endif
