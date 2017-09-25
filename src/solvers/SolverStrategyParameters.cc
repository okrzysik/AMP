#include "SolverStrategyParameters.h"

#include <utility>


namespace AMP {
namespace Solver {

SolverStrategyParameters::SolverStrategyParameters() { d_name = "SolverStrategyParameters"; }

SolverStrategyParameters::SolverStrategyParameters( AMP::shared_ptr<AMP::Database> db )
    : d_db( std::move( db ) )
{
    d_name = "SolverStrategyParameters";
}

SolverStrategyParameters::~SolverStrategyParameters() = default;
} // namespace Solver
} // namespace AMP
