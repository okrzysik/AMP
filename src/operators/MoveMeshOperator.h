
#ifndef included_AMP_MoveMeshOperator
#define included_AMP_MoveMeshOperator

#include "AMP/operators/Operator.h"

namespace AMP {
namespace Operator {

class MoveMeshOperator : public Operator
{
public:
    explicit MoveMeshOperator( std::shared_ptr<const OperatorParameters> params );

    virtual ~MoveMeshOperator() {}

    std::string type() const override { return "MoveMeshOperator"; }

    void setVariable( AMP::LinearAlgebra::Variable::shared_ptr var );

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override;

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

protected:
    AMP::LinearAlgebra::Variable::shared_ptr d_var;
    AMP::LinearAlgebra::Vector::shared_ptr d_prevDisp;
};
} // namespace Operator
} // namespace AMP

#endif
