
#ifndef included_AMP_OnePointOperator
#define included_AMP_OnePointOperator

#include "AMP/operators/Operator.h"

namespace AMP::Operator {

class OnePointOperator : public Operator
{
public:
    explicit OnePointOperator( std::shared_ptr<const OperatorParameters> params )
        : Operator( params )
    {
        d_constant = 0.0;
    }

    std::string type() const override { return "OnePointOperator"; }

    double getConstant() { return d_constant; }

protected:
    double d_constant;

private:
};
} // namespace AMP::Operator

#endif
