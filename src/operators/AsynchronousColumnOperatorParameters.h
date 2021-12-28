#ifndef included_AMP_AsynchronousColumnOperatorParameters
#define included_AMP_AsynchronousColumnOperatorParameters

#include "AMP/operators/OperatorParameters.h"

namespace AMP::Operator {

class AsynchronousColumnOperatorParameters : public OperatorParameters
{
public:
    explicit AsynchronousColumnOperatorParameters( std::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
    }
};
} // namespace AMP::Operator

#endif
