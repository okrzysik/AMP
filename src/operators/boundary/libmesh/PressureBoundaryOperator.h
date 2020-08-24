
#ifndef included_AMP_PressureBoundaryOperator
#define included_AMP_PressureBoundaryOperator

#include "AMP/operators/boundary/libmesh/TractionBoundaryOperator.h"

namespace AMP {
namespace Operator {

class PressureBoundaryOperator : public BoundaryOperator
{
public:
    explicit PressureBoundaryOperator( const std::shared_ptr<OperatorParameters> &params );

    virtual ~PressureBoundaryOperator() {}

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr,
                AMP::LinearAlgebra::Vector::shared_ptr r ) override
    {
        AMP::LinearAlgebra::Vector::shared_ptr nullVec;
        d_tractionOp->apply( nullVec, r );
        r->scale( -1.0, *r );
    }

    void addRHScorrection( AMP::LinearAlgebra::Vector::shared_ptr rhs ) override
    {
        d_tractionOp->addRHScorrection( rhs );
    }

protected:
    std::shared_ptr<AMP::Operator::TractionBoundaryOperator> d_tractionOp;
};
} // namespace Operator
} // namespace AMP

#endif
