#include "AMP/solvers/KrylovSolverParameters.h"

namespace AMP::Solver {

KrylovSolverParameters::KrylovSolverParameters( std::shared_ptr<AMP::Database> db )
    : SolverStrategyParameters( db )
{
}
} // namespace AMP::Solver
