#ifndef included_AMP_CoupledFlowFrapconOperatorParameters
#define included_AMP_CoupledFlowFrapconOperatorParameters

#include "AMP/operators/ColumnOperatorParameters.h"
#include "AMP/operators/Operator.h"
#include <memory>

#include <vector>

namespace AMP::Operator {

/**
  A class that encapsulates the parameters required to construct
  the composite Operator operator.
  @see ColumnOperator
  */
class CoupledFlowFrapconOperatorParameters : public ColumnOperatorParameters
{
public:
    explicit CoupledFlowFrapconOperatorParameters( std::shared_ptr<AMP::Database> db )
        : ColumnOperatorParameters( db )
    {
    }

    virtual ~CoupledFlowFrapconOperatorParameters() {}

    std::shared_ptr<Operator> d_Map3to1;

    std::shared_ptr<Operator> d_Map1to3;

    std::shared_ptr<Operator> d_FlowOperator;
};
} // namespace AMP::Operator


#endif
