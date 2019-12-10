
#ifndef included_AMP_CoupledChannelToCladMapOperatorParameters
#define included_AMP_CoupledChannelToCladMapOperatorParameters

/*AMP files */
#include "AMP/operators/Operator.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/operators/subchannel/SubchannelPhysicsModel.h"
/*Boost files */
#include "AMP/utils/shared_ptr.h"

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
    explicit CoupledChannelToCladMapOperatorParameters( AMP::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
    }

    virtual ~CoupledChannelToCladMapOperatorParameters() {}

    AMP::shared_ptr<AMP::Operator::Operator> d_thermalMapOperator;
    AMP::shared_ptr<AMP::Operator::Operator> d_densityMapOperator;
    AMP::shared_ptr<AMP::Operator::Operator> d_mapOperator;
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_variable;
    AMP::shared_ptr<AMP::Mesh::Mesh> d_subchannelMesh;
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_vector;
    AMP::shared_ptr<AMP::Operator::SubchannelPhysicsModel> d_subchannelPhysicsModel;
};
} // namespace Operator
} // namespace AMP


#endif
