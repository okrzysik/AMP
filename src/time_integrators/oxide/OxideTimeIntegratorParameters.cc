#include "time_integrators/oxide/OxideTimeIntegratorParameters.h"

namespace AMP{
namespace TimeIntegrator{


/************************************************************************
* Constructor and destructor for TimeIntegrator.                        *
************************************************************************/
OxideTimeIntegratorParameters::OxideTimeIntegratorParameters(const boost::shared_ptr<AMP::Database> db) :
    TimeIntegratorParameters( db )
{
    depth = -1.0;
}



}
}

