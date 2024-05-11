#include "AMP/time_integrators/TimeIntegratorParameters.h"

namespace AMP::TimeIntegrator {

TimeIntegratorParameters::TimeIntegratorParameters( std::shared_ptr<AMP::Database> db )
    : Operator::OperatorParameters( db )
{
}

TimeIntegratorParameters::~TimeIntegratorParameters() = default;
} // namespace AMP::TimeIntegrator
