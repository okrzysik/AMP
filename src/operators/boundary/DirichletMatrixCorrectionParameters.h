
#ifndef included_AMP_DirichletMatrixCorrectionParameters
#define included_AMP_DirichletMatrixCorrectionParameters

#include "AMP/operators/boundary/LinearBoundaryOperatorParameters.h"

namespace AMP {
namespace Operator {

class DirichletMatrixCorrectionParameters : public LinearBoundaryOperatorParameters
{
public:
    explicit DirichletMatrixCorrectionParameters( std::shared_ptr<AMP::Database> db )
        : LinearBoundaryOperatorParameters( db )
    {
    }

    virtual ~DirichletMatrixCorrectionParameters() {}

    // This must be a simple variable not a dual or multivariable
    AMP::LinearAlgebra::Variable::shared_ptr d_variable;
};
} // namespace Operator
} // namespace AMP

#endif
