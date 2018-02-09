
#ifndef included_AMP_RobinVectorCorrectionParameters
#define included_AMP_RobinVectorCorrectionParameters

#include "AMP/operators/OperatorParameters.h"

namespace AMP {
namespace Operator {

class RobinVectorCorrectionParameters : public OperatorParameters
{
public:
    explicit RobinVectorCorrectionParameters( const AMP::shared_ptr<AMP::Database> &db )
        : OperatorParameters( db )
    {
    }

    virtual ~RobinVectorCorrectionParameters() {}

    AMP::LinearAlgebra::Variable::shared_ptr d_variable;
};
} // namespace Operator
} // namespace AMP

#endif
