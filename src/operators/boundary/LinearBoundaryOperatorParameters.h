
#ifndef included_AMP_LinearBoundaryOperatorParameters
#define included_AMP_LinearBoundaryOperatorParameters

#include "AMP/matrices/Matrix.h"
#include "AMP/operators/OperatorParameters.h"

namespace AMP {
namespace Operator {

class LinearBoundaryOperatorParameters : public OperatorParameters
{
public:
    explicit LinearBoundaryOperatorParameters( AMP::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
    }

    virtual ~LinearBoundaryOperatorParameters() {}

    AMP::LinearAlgebra::Matrix::shared_ptr d_inputMatrix;
};
} // namespace Operator
} // namespace AMP

#endif
