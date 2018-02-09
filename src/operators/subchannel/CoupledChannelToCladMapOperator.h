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
        const AMP::shared_ptr<CoupledChannelToCladMapOperatorParameters> &params );

    virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override
    {
        return d_flowVariable;
    }

    virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override
    {
        return d_thermalMapOperator->getOutputVariable();
    }

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr r ) override;

    virtual ~CoupledChannelToCladMapOperator() {}

protected:
private:
    AMP::LinearAlgebra::Variable::shared_ptr d_flowVariable;

    AMP::LinearAlgebra::Vector::shared_ptr d_subchannelTemperature;
    AMP::LinearAlgebra::Vector::shared_ptr d_subchannelDensity;

    AMP::shared_ptr<AMP::Operator::Operator> d_thermalMapOperator;
    AMP::shared_ptr<AMP::Operator::Operator> d_densityMapOperator;

    AMP::shared_ptr<SubchannelPhysicsModel> d_subchannelPhysicsModel;
};
} // namespace Operator
} // namespace AMP

#endif
