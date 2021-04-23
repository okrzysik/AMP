#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"

namespace AMP {
namespace Solver {

PetscKrylovSolverParameters::PetscKrylovSolverParameters( std::shared_ptr<AMP::Database> db )
    : SolverStrategyParameters( db )
{
}
} // namespace Solver
} // namespace AMP
