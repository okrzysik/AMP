#ifndef included_AMP_LinearCoupledFlowOperator
#define included_AMP_LinearCoupledFlowOperator

#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Vector.h"
#include "LinearCoupledFlowOperatorParameters.h"
#include <vector>

namespace AMP {
namespace Operator {

class LinearCoupledFlowOperator : public Operator
{
public:
    explicit LinearCoupledFlowOperator( std::shared_ptr<const OperatorParameters> params )
        : Operator( params )
    {
        (void) params;
    }

    virtual ~LinearCoupledFlowOperator() {}

    std::string type() const override { return "LinearCoupledFlowOperator"; }

    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                        AMP::LinearAlgebra::Vector::const_shared_ptr u,
                        AMP::LinearAlgebra::Vector::shared_ptr r,
                        const double a = -1.0,
                        const double b = 1.0 );

    virtual void reset( std::shared_ptr<const OperatorParameters> params );

    virtual void append( std::shared_ptr<Operator> op );

    virtual std::shared_ptr<AMP::LinearAlgebra::Variable> getOutputVariable();

    virtual std::shared_ptr<AMP::LinearAlgebra::Variable> getInputVariable( int varId = -1 );

protected:
    std::vector<std::shared_ptr<Operator>> d_Operators;

private:
};
} // namespace Operator
} // namespace AMP

#endif
