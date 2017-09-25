#ifndef included_AMP_LinearCoupledFlowOperator
#define included_AMP_LinearCoupledFlowOperator

#include "LinearCoupledFlowOperatorParameters.h"
#include "utils/Utilities.h"
#include "vectors/Vector.h"
#include <vector>

namespace AMP {
namespace Operator {

class LinearCoupledFlowOperator : public Operator
{
public:
    explicit LinearCoupledFlowOperator( const AMP::shared_ptr<OperatorParameters> &params )
        : Operator( params )
    {
        (void) params;
    }

    virtual ~LinearCoupledFlowOperator() {}

    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                        AMP::LinearAlgebra::Vector::const_shared_ptr u,
                        AMP::LinearAlgebra::Vector::shared_ptr r,
                        const double a = -1.0,
                        const double b = 1.0 );

    virtual void reset( const AMP::shared_ptr<OperatorParameters> &params );

    virtual void append( AMP::shared_ptr<Operator> op );

    virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable();

    virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable( int varId = -1 );

protected:
    std::vector<AMP::shared_ptr<Operator>> d_Operators;

private:
};
} // namespace Operator
} // namespace AMP

#endif
