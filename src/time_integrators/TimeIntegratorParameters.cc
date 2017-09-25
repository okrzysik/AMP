#include "TimeIntegratorParameters.h"

namespace AMP {
namespace TimeIntegrator {

TimeIntegratorParameters::TimeIntegratorParameters( const AMP::shared_ptr<AMP::Database> db )
{
    d_db = db;
    d_ic_vector.reset();
    d_operator.reset();
}

TimeIntegratorParameters::~TimeIntegratorParameters() = default;
} // namespace TimeIntegrator
} // namespace AMP
