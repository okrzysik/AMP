#ifndef included_AMP_AsynchronousColumnOperatorParameters
#define included_AMP_AsynchronousColumnOperatorParameters

#include "AMP/operators/OperatorParameters.h"

namespace AMP {
namespace Operator {

class AsynchronousColumnOperatorParameters : public OperatorParameters
{
public:
    explicit AsynchronousColumnOperatorParameters( AMP::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
    }
};
} // namespace Operator
} // namespace AMP

#endif
