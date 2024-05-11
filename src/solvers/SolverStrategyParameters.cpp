#include "AMP/solvers/SolverStrategyParameters.h"

#include <utility>


namespace AMP::Solver {

SolverStrategyParameters::SolverStrategyParameters() { d_name = "SolverStrategyParameters"; }

SolverStrategyParameters::SolverStrategyParameters( std::shared_ptr<AMP::Database> db )
    : ParameterBase( std::move( db ) )
{
    d_name = "SolverStrategyParameters";
}

SolverStrategyParameters::~SolverStrategyParameters() = default;
} // namespace AMP::Solver
