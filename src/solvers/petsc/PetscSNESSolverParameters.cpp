#include "AMP/solvers/petsc/PetscSNESSolverParameters.h"

namespace AMP::Solver {

PetscSNESSolverParameters::PetscSNESSolverParameters( std::shared_ptr<AMP::Database> db )
    : SolverStrategyParameters( db )
{
}
} // namespace AMP::Solver
