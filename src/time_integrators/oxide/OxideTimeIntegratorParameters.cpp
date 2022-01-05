#include "AMP/time_integrators/oxide/OxideTimeIntegratorParameters.h"

namespace AMP::TimeIntegrator {


/************************************************************************
 * Constructor and destructor for TimeIntegrator.                        *
 ************************************************************************/
OxideTimeIntegratorParameters::OxideTimeIntegratorParameters( std::shared_ptr<AMP::Database> db )
    : TimeIntegratorParameters( db )
{
    depth = -1.0;
}
} // namespace AMP::TimeIntegrator
