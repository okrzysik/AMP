#ifndef included_AMP_CoupledChannelToCladMapOperator
#define included_AMP_CoupledChannelToCladMapOperator

#include "AMP/operators/Operator.h"
#include "AMP/operators/subchannel/CoupledChannelToCladMapOperatorParameters.h"
#include "AMP/operators/subchannel/SubchannelPhysicsModel.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Vector.h"
#include <vector>

namespace AMP {
namespace Operator {

class CoupledChannelToCladMapOperator : public Operator
{
public:
    explicit CoupledChannelToCladMapOperator(
        const std::shared_ptr<CoupledChannelToCladMapOperatorParameters> &params );

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override { return d_flowVariable; }

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override
    {
        return d_thermalMapOperator->getOutputVariable();
    }

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr r ) override;

    virtual ~CoupledChannelToCladMapOperator() {}

    std::string type() const override { return "CoupledChannelToCladMapOperator"; }

protected:
private:
    AMP::LinearAlgebra::Variable::shared_ptr d_flowVariable;

    AMP::LinearAlgebra::Vector::shared_ptr d_subchannelTemperature;
    AMP::LinearAlgebra::Vector::shared_ptr d_subchannelDensity;

    std::shared_ptr<AMP::Operator::Operator> d_thermalMapOperator;
    std::shared_ptr<AMP::Operator::Operator> d_densityMapOperator;

    std::shared_ptr<SubchannelPhysicsModel> d_subchannelPhysicsModel;
};
} // namespace Operator
} // namespace AMP

#endif
