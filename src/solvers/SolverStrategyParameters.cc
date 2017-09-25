#include "SolverStrategyParameters.h"


namespace AMP {
namespace Solver {

SolverStrategyParameters::SolverStrategyParameters() { d_name = "SolverStrategyParameters"; }

SolverStrategyParameters::SolverStrategyParameters( const AMP::shared_ptr<AMP::Database> &db )
    : d_db( db )
{
    d_name = "SolverStrategyParameters";
}

SolverStrategyParameters::~SolverStrategyParameters() {}
} // namespace Solver
} // namespace AMP
