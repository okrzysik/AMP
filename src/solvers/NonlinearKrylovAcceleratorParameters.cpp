
#include "NonlinearKrylovAcceleratorParameters.h"

namespace AMP::Solver {

NonlinearKrylovAcceleratorParameters::NonlinearKrylovAcceleratorParameters() = default;

NonlinearKrylovAcceleratorParameters::NonlinearKrylovAcceleratorParameters(
    std::shared_ptr<AMP::Database> database )
    : SolverStrategyParameters( database )
{
}

NonlinearKrylovAcceleratorParameters::~NonlinearKrylovAcceleratorParameters() = default;
} // namespace AMP::Solver
