#include "solvers/petsc/PetscKrylovSolverParameters.h"

namespace AMP {
namespace Solver {

PetscKrylovSolverParameters::PetscKrylovSolverParameters( const AMP::shared_ptr<AMP::Database> db )
    : SolverStrategyParameters( db )
{
}
} // namespace Solver
} // namespace AMP
