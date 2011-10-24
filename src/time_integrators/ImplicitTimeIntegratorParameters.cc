#include "ImplicitTimeIntegratorParameters.h"

namespace AMP{
namespace TimeIntegrator{

ImplicitTimeIntegratorParameters::ImplicitTimeIntegratorParameters( boost::shared_ptr<AMP::Database> db)
:TimeIntegratorParameters(db)
{
  d_solver.reset();
}

ImplicitTimeIntegratorParameters::~ImplicitTimeIntegratorParameters()
{
}

}
}

