
#ifndef included_AMP_ColumnOperatorParameters
#define included_AMP_ColumnOperatorParameters

/*AMP files */
#include "operators/OperatorParameters.h"

/*Boost files */
#include "utils/shared_ptr.h"

#include <vector>

namespace AMP {
namespace Operator {

/**
  A class that encapsulates the parameters required to construct
  the composite Operator operator.
  @see ColumnOperator
  */
class ColumnOperatorParameters : public OperatorParameters {
public:
    explicit ColumnOperatorParameters( const AMP::shared_ptr<AMP::Database> &db )
        : OperatorParameters( db )
    {
    }

    virtual ~ColumnOperatorParameters() {}

    std::vector<AMP::shared_ptr<OperatorParameters>> d_OperatorParameters;
};
}
}


#endif
