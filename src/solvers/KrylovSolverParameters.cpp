#include "AMP/solvers/KrylovSolverParameters.h"

namespace AMP {
namespace Solver {

KrylovSolverParameters::KrylovSolverParameters( const std::shared_ptr<AMP::Database> db )
    : SolverStrategyParameters( db )
{
}
} // namespace Solver
} // namespace AMP
