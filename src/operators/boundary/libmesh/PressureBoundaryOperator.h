
#ifndef included_AMP_PressureBoundaryOperator
#define included_AMP_PressureBoundaryOperator

#include "operators/boundary/libmesh/TractionBoundaryOperator.h"

namespace AMP {
namespace Operator {

class PressureBoundaryOperator : public BoundaryOperator {
public:
    explicit PressureBoundaryOperator( const AMP::shared_ptr<OperatorParameters> &params );

    virtual ~PressureBoundaryOperator() {}

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr,
                AMP::LinearAlgebra::Vector::shared_ptr r ) override
    {
        AMP::LinearAlgebra::Vector::shared_ptr nullVec;
        d_tractionOp->apply( nullVec, r );
        r->scale( -1.0 );
    }

    void addRHScorrection( AMP::LinearAlgebra::Vector::shared_ptr rhs )
    {
        d_tractionOp->addRHScorrection( rhs );
    }

protected:
    AMP::shared_ptr<AMP::Operator::TractionBoundaryOperator> d_tractionOp;
};
}
}

#endif
