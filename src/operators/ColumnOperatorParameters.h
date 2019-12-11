#ifndef included_AMP_ColumnOperatorParameters
#define included_AMP_ColumnOperatorParameters

#include "AMP/operators/OperatorParameters.h"
#include "AMP/utils/shared_ptr.h"

#include <vector>

namespace AMP {
namespace Operator {

/**
  A class that encapsulates the parameters required to construct
  the composite Operator operator.
  @see ColumnOperator
  */
class ColumnOperatorParameters : public OperatorParameters
{
public:
    explicit ColumnOperatorParameters( AMP::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
    }

    virtual ~ColumnOperatorParameters() {}

    std::vector<AMP::shared_ptr<OperatorParameters>> d_OperatorParameters;
};
} // namespace Operator
} // namespace AMP


#endif
