
#ifndef included_AMP_CoupledFlowFrapconOperatorParameters
#define included_AMP_CoupledFlowFrapconOperatorParameters

/*AMP files */
#include "AMP/operators/ColumnOperatorParameters.h"
#include "AMP/operators/Operator.h"
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
class CoupledFlowFrapconOperatorParameters : public ColumnOperatorParameters
{
public:
    explicit CoupledFlowFrapconOperatorParameters( const AMP::shared_ptr<AMP::Database> &db )
        : ColumnOperatorParameters( db )
    {
    }

    virtual ~CoupledFlowFrapconOperatorParameters() {}

    AMP::shared_ptr<Operator> d_Map3to1;

    AMP::shared_ptr<Operator> d_Map1to3;

    AMP::shared_ptr<Operator> d_FlowOperator;
};
} // namespace Operator
} // namespace AMP


#endif
