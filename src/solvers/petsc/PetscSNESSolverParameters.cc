#include "AMP/solvers/petsc/PetscSNESSolverParameters.h"

namespace AMP {
namespace Solver {

PetscSNESSolverParameters::PetscSNESSolverParameters( AMP::shared_ptr<AMP::Database> db )
    : SolverStrategyParameters( db )
{
}
} // namespace Solver
} // namespace AMP
