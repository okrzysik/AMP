#include "time_integrators/sundials/IDATimeIntegratorParameters.h"

namespace AMP {
namespace TimeIntegrator {

IDATimeIntegratorParameters::IDATimeIntegratorParameters( const AMP::shared_ptr<AMP::Database> db )
    : TimeIntegratorParameters( db )

{
    // Question: does this automatically call the constructor of its parent class?
    /*
    d_db = db;
    d_ic_vector.reset();:
    d_operator.reset();
     */
    d_ic_vector_prime.reset();
}

IDATimeIntegratorParameters::~IDATimeIntegratorParameters() = default;
} // namespace TimeIntegrator
} // namespace AMP
