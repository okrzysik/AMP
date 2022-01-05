#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"

namespace AMP::Solver {

PetscKrylovSolverParameters::PetscKrylovSolverParameters( std::shared_ptr<AMP::Database> db )
    : SolverStrategyParameters( db )
{
}
} // namespace AMP::Solver
