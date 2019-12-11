
#ifndef included_AMP_FlowFrapconJacobianParameters
#define included_AMP_FlowFrapconJacobianParameters

#include "AMP/operators/OperatorParameters.h"

namespace AMP {
namespace Operator {

class FlowFrapconJacobianParameters : public OperatorParameters
{
public:
    explicit FlowFrapconJacobianParameters( std::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
    }

    virtual ~FlowFrapconJacobianParameters() {}

    AMP::LinearAlgebra::Variable::shared_ptr d_variable;

    AMP::LinearAlgebra::Vector::shared_ptr d_frozenSolution;
};
} // namespace Operator
} // namespace AMP

#endif
