#ifndef included_AMP_CoupledChannelToCladMapOperatorParameters
#define included_AMP_CoupledChannelToCladMapOperatorParameters

#include "AMP/operators/Operator.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/operators/subchannel/SubchannelPhysicsModel.h"
#include <memory>

#include <vector>

namespace AMP {
namespace Operator {

/**
  A class that encapsulates the parameters required to construct
  the composite Operator operator.
  @see ColumnOperator
  */
class CoupledChannelToCladMapOperatorParameters : public OperatorParameters
{
public:
    explicit CoupledChannelToCladMapOperatorParameters( std::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
    }

    virtual ~CoupledChannelToCladMapOperatorParameters() {}

    std::shared_ptr<AMP::Operator::Operator> d_thermalMapOperator;
    std::shared_ptr<AMP::Operator::Operator> d_densityMapOperator;
    std::shared_ptr<AMP::Operator::Operator> d_mapOperator;
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_variable;
    std::shared_ptr<AMP::Mesh::Mesh> d_subchannelMesh;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_vector;
    std::shared_ptr<AMP::Operator::SubchannelPhysicsModel> d_subchannelPhysicsModel;
};
} // namespace Operator
} // namespace AMP


#endif
