#ifndef included_AMP_ColumnOperatorParameters
#define included_AMP_ColumnOperatorParameters

#include "AMP/operators/OperatorParameters.h"
#include <memory>

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
    explicit ColumnOperatorParameters( std::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
    }

    virtual ~ColumnOperatorParameters() {}

    std::vector<std::shared_ptr<OperatorParameters>> d_OperatorParameters;
};
} // namespace Operator
} // namespace AMP


#endif
