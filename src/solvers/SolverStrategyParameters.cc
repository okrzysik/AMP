#include "SolverStrategyParameters.h"


namespace AMP {
namespace Solver {

SolverStrategyParameters::SolverStrategyParameters()
{
   d_name="SolverStrategyParameters";
}

  SolverStrategyParameters::SolverStrategyParameters(const boost::shared_ptr<AMP::Database> &db)
:d_db(db)
{
   d_name     = "SolverStrategyParameters";
}

SolverStrategyParameters::~SolverStrategyParameters()
{
}

}
}

