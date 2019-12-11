#include "TimeIntegratorParameters.h"

namespace AMP {
namespace TimeIntegrator {

TimeIntegratorParameters::TimeIntegratorParameters( const std::shared_ptr<AMP::Database> db )
    : d_db( db )
{
}

TimeIntegratorParameters::~TimeIntegratorParameters() = default;
} // namespace TimeIntegrator
} // namespace AMP
