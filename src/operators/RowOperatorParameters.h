#ifndef included_AMP_RowOperatorParameters
#define included_AMP_RowOperatorParameters

#include "AMP/operators/Operator.h"
#include "AMP/operators/OperatorParameters.h"
#include <memory>

namespace AMP {
namespace Operator {

class RowOperatorParameters : public OperatorParameters
{
public:
    explicit RowOperatorParameters( std::shared_ptr<AMP::Database> db ) : OperatorParameters( db )
    {
    }

    virtual ~RowOperatorParameters(){};


    std::vector<std::shared_ptr<AMP::Operator>> d_Operator;

    std::vector<std::shared_ptr<AMP::OperatorParameters>> d_OperatorParameters;

    std::vector<double> scalea;
};
} // namespace Operator
} // namespace AMP

#endif
