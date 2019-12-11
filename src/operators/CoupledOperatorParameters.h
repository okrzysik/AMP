#ifndef included_AMP_CoupledOperatorParameters
#define included_AMP_CoupledOperatorParameters

#include "AMP/operators/ColumnOperatorParameters.h"
#include "AMP/operators/Operator.h"
#include <memory>

#include <vector>

namespace AMP {
namespace Operator {

/**
  A class that encapsulates the parameters required to construct
  the composite Operator operator.
  @see ColumnOperator
  */
class CoupledOperatorParameters : public ColumnOperatorParameters
{
public:
    explicit CoupledOperatorParameters( std::shared_ptr<AMP::Database> db )
        : ColumnOperatorParameters( db )
    {
    }

    virtual ~CoupledOperatorParameters() {}

    std::shared_ptr<Operator> d_NodeToGaussPointOperator;

    std::shared_ptr<Operator> d_CopyOperator;

    std::shared_ptr<Operator> d_MapOperator;

    std::shared_ptr<Operator> d_BVPOperator;
};
} // namespace Operator
} // namespace AMP


#endif
