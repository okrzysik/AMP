#include "ImplicitTimeIntegratorParameters.h"

namespace AMP {
namespace TimeIntegrator {

ImplicitTimeIntegratorParameters::ImplicitTimeIntegratorParameters(
    std::shared_ptr<AMP::Database> db )
    : TimeIntegratorParameters( db )
{
    d_solver.reset();
}

ImplicitTimeIntegratorParameters::~ImplicitTimeIntegratorParameters() = default;
} // namespace TimeIntegrator
} // namespace AMP
