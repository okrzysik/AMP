#include "TimeIntegratorParameters.h"

namespace AMP {
namespace TimeIntegrator {

TimeIntegratorParameters::TimeIntegratorParameters( std::shared_ptr<AMP::Database> db ) : d_db( db )
{
}

TimeIntegratorParameters::~TimeIntegratorParameters() = default;
} // namespace TimeIntegrator
} // namespace AMP
