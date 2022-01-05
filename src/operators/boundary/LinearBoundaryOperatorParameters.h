
#ifndef included_AMP_LinearBoundaryOperatorParameters
#define included_AMP_LinearBoundaryOperatorParameters

#include "AMP/matrices/Matrix.h"
#include "AMP/operators/OperatorParameters.h"

namespace AMP::Operator {

class LinearBoundaryOperatorParameters : public OperatorParameters
{
public:
    explicit LinearBoundaryOperatorParameters( std::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
    }

    virtual ~LinearBoundaryOperatorParameters() {}

    std::shared_ptr<AMP::LinearAlgebra::Matrix> d_inputMatrix;
};
} // namespace AMP::Operator

#endif
