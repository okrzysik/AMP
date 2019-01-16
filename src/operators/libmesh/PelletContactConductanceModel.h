#ifndef included_AMP_PelletContactConductanceModel
#define included_AMP_PelletContactConductanceModel

#include <string>
#include <vector>

#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/boundary/libmesh/RobinPhysicsModel.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/utils/shared_ptr.h"

namespace AMP {
namespace Operator {

class PelletContactConductanceModel : public RobinPhysicsModel
{
public:
    explicit PelletContactConductanceModel(
        const AMP::shared_ptr<RobinPhysicsModelParameters> &params );

    virtual ~PelletContactConductanceModel() {}

    void getConductance( std::vector<double> &beta,
                         std::vector<double> &gamma,
                         const std::vector<std::vector<double>> &arguments ) override;

protected:
    unsigned int d_nTransportModels; /**< Number of Transport Models. */
    std::vector<AMP::shared_ptr<DiffusionTransportModel>> d_transportModels;

private:
};
} // namespace Operator
} // namespace AMP

#endif
