
#ifndef included_AMP_FlowFrapconOperatorParameters
#define included_AMP_FlowFrapconOperatorParameters

#include "AMP/operators/OperatorParameters.h"
#include "AMP/operators/subchannel/FlowFrapconJacobianParameters.h"

namespace AMP {
namespace Operator {

class FlowFrapconOperatorParameters : public OperatorParameters
{
public:
    explicit FlowFrapconOperatorParameters( std::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
    }

    virtual ~FlowFrapconOperatorParameters() {}

    std::shared_ptr<AMP::LinearAlgebra::Variable> d_variable;
};
} // namespace Operator
} // namespace AMP

#endif
