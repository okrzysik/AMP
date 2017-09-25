
#include "NonlinearKrylovAcceleratorParameters.h"

namespace AMP {
namespace Solver {

NonlinearKrylovAcceleratorParameters::NonlinearKrylovAcceleratorParameters() = default;

NonlinearKrylovAcceleratorParameters::NonlinearKrylovAcceleratorParameters(
    const AMP::shared_ptr<AMP::Database> &database )
    : SolverStrategyParameters( database )
{
}

NonlinearKrylovAcceleratorParameters::~NonlinearKrylovAcceleratorParameters() = default;
} // namespace Solver
} // namespace AMP
