#ifndef included_AMP_ConvectiveHeatCoefficient
#define included_AMP_ConvectiveHeatCoefficient

#include <iostream>
#include <string>
#include <vector>

#include "operators/ElementPhysicsModelFactory.h"
#include "operators/boundary/libmesh/RobinPhysicsModel.h"
#include "utils/shared_ptr.h"

namespace AMP {
namespace Operator {

class ConvectiveHeatCoefficient : public RobinPhysicsModel
{
public:
    explicit ConvectiveHeatCoefficient(
        const AMP::shared_ptr<RobinPhysicsModelParameters> &params );

    virtual ~ConvectiveHeatCoefficient() {}

    void getConductance( std::vector<double> &beta,
                         std::vector<double> &gamma,
                         const std::vector<std::vector<double>> &arguments );

    AMP::Materials::Material::shared_ptr getMaterial() { return d_material; }
    AMP::shared_ptr<AMP::Materials::Property<double>> getProperty() { return d_property; }

protected:
    AMP::Materials::Material::shared_ptr d_material;

    AMP::shared_ptr<AMP::Materials::Property<double>> d_property;

    std::vector<std::string> d_argNames;

    std::vector<double> d_defaults;

private:
};
} // namespace Operator
} // namespace AMP

#endif
