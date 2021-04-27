#ifndef included_AMP_PetscKrylovSolverParameters
#define included_AMP_PetscKrylovSolverParameters

#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/utils/Database.h"
#include <memory>


namespace AMP {
namespace Solver {


/**
 * Class PetscKrylovSolverParameters provides a uniform mechanism to pass
 * initialization parameters to the PetscKrylovSolver solver. It contains
 * shared pointers to a database object and a preconditioner. All member variables are public.
 */
class PetscKrylovSolverParameters : public SolverStrategyParameters
{
public:
    PetscKrylovSolverParameters() {}
    explicit PetscKrylovSolverParameters( std::shared_ptr<AMP::Database> db );
    virtual ~PetscKrylovSolverParameters() {}

    AMP_MPI d_comm;

    std::shared_ptr<AMP::Solver::SolverStrategy> d_pPreconditioner;

protected:
private:
};


} // namespace Solver
} // namespace AMP

#endif
