
#ifndef included_LinearBoundaryOperatorParameters
#define included_LinearBoundaryOperatorParameters

#include "matrices/Matrix.h"
#include "operators/OperatorParameters.h"

namespace AMP {
namespace Operator {

class LinearBoundaryOperatorParameters : public OperatorParameters
{
public:
    explicit LinearBoundaryOperatorParameters( const AMP::shared_ptr<AMP::Database> &db )
        : OperatorParameters( db )
    {
    }

    virtual ~LinearBoundaryOperatorParameters() {}

    AMP::LinearAlgebra::Matrix::shared_ptr d_inputMatrix;
};
} // namespace Operator
} // namespace AMP

#endif
