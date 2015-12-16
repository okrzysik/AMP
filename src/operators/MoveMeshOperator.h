
#ifndef included_AMP_MoveMeshOperator
#define included_AMP_MoveMeshOperator

#include "operators/Operator.h"

namespace AMP {
namespace Operator {

class MoveMeshOperator : public Operator
{
public:
    explicit MoveMeshOperator( const AMP::shared_ptr<OperatorParameters> &params );

    virtual ~MoveMeshOperator() {}

    void setVariable( AMP::LinearAlgebra::Variable::shared_ptr var );

    virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override;

    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                        AMP::LinearAlgebra::Vector::shared_ptr f ) override;

protected:
    AMP::LinearAlgebra::Variable::shared_ptr d_var;
    AMP::LinearAlgebra::Vector::shared_ptr d_prevDisp;
};
}
}

#endif
