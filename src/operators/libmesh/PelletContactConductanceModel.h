#ifndef included_AMP_PelletContactConductanceModel
#define included_AMP_PelletContactConductanceModel

#include <string>
#include <vector>

#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/boundary/libmesh/RobinPhysicsModel.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include <memory>

namespace AMP::Operator {

class PelletContactConductanceModel : public RobinPhysicsModel
{
public:
    explicit PelletContactConductanceModel(
        std::shared_ptr<const RobinPhysicsModelParameters> params );

    virtual ~PelletContactConductanceModel() {}

    void getConductance( std::vector<double> &beta,
                         std::vector<double> &gamma,
                         const std::vector<std::vector<double>> &arguments ) override;

protected:
    unsigned int d_nTransportModels; /**< Number of Transport Models. */
    std::vector<std::shared_ptr<DiffusionTransportModel>> d_transportModels;

private:
};
} // namespace AMP::Operator

#endif
