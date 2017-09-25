#ifndef included_OxideTimeIntegratorParameters
#define included_OxideTimeIntegratorParameters

#include "ampmesh/Mesh.h"
#include "time_integrators/TimeIntegratorParameters.h"

namespace AMP {
namespace TimeIntegrator {


/*!
  @brief This class contains the parmeters to create the problem for the oxide growth
 */
class OxideTimeIntegratorParameters : public TimeIntegratorParameters
{
public:
    explicit OxideTimeIntegratorParameters( const AMP::shared_ptr<AMP::Database> db );

    // Surface mesh for the oxide calculation
    AMP::Mesh::Mesh::shared_ptr d_mesh;

    // Temperature vector (K)
    AMP::LinearAlgebra::Vector::shared_ptr d_temp;

    // Surface thickness (m)
    double depth;
};
} // namespace TimeIntegrator
} // namespace AMP

#endif
