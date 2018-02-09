#ifndef included_RowOperatorParameters
#define included_RowOperatorParameters

#include "AMP/operators/Operator.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/utils/shared_ptr.h"

namespace AMP {
namespace Operator {

class RowOperatorParameters : public OperatorParameters
{
public:
    explicit RowOperatorParameters( const AMP::shared_ptr<AMP::Database> &db )
        : OperatorParameters( db )
    {
    }

    virtual ~RowOperatorParameters(){};


    std::vector<AMP::shared_ptr<AMP::Operator>> d_Operator;

    std::vector<AMP::shared_ptr<AMP::OperatorParameters>> d_OperatorParameters;

    std::vector<double> scalea;
};
} // namespace Operator
} // namespace AMP

#endif
