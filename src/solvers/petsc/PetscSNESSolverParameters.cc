#include "AMP/solvers/petsc/PetscSNESSolverParameters.h"

namespace AMP {
namespace Solver {

PetscSNESSolverParameters::PetscSNESSolverParameters( const AMP::shared_ptr<AMP::Database> &db )
    : SolverStrategyParameters( db )
{
}
} // namespace Solver
} // namespace AMP
