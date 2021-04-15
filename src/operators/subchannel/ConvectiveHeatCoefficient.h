#ifndef included_AMP_ConvectiveHeatCoefficient
#define included_AMP_ConvectiveHeatCoefficient

#include <string>
#include <vector>

#include "AMP/materials/Material.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/boundary/libmesh/RobinPhysicsModel.h"
#include <memory>


namespace AMP {
namespace Operator {

class ConvectiveHeatCoefficient : public RobinPhysicsModel
{
public:
    explicit ConvectiveHeatCoefficient(
        const std::shared_ptr<RobinPhysicsModelParameters> &params );

    virtual ~ConvectiveHeatCoefficient() {}

    void getConductance( std::vector<double> &beta,
                         std::vector<double> &gamma,
                         const std::vector<std::vector<double>> &arguments ) override;

    std::shared_ptr<AMP::Materials::Material> getMaterial() { return d_material; }
    std::shared_ptr<AMP::Materials::Property<double>> getProperty() { return d_property; }

protected:
    std::shared_ptr<AMP::Materials::Material> d_material;

    std::shared_ptr<AMP::Materials::Property<double>> d_property;

    std::vector<std::string> d_argNames;

    std::vector<double> d_defaults;

private:
};
} // namespace Operator
} // namespace AMP

#endif
