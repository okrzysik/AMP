
#include "NonlinearKrylovAcceleratorParameters.h"

namespace AMP {
namespace Solver {

NonlinearKrylovAcceleratorParameters :: NonlinearKrylovAcceleratorParameters()
{
}

NonlinearKrylovAcceleratorParameters :: NonlinearKrylovAcceleratorParameters(const boost::shared_ptr<AMP::Database> &database)
: SolverStrategyParameters(database)
{
}

NonlinearKrylovAcceleratorParameters :: ~NonlinearKrylovAcceleratorParameters()
{
}
  
}
}

