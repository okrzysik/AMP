
#ifndef included_AMP_DirichletVectorCorrectionParameters
#define included_AMP_DirichletVectorCorrectionParameters

#include "AMP/operators/OperatorParameters.h"

namespace AMP {
namespace Operator {

class DirichletVectorCorrectionParameters : public OperatorParameters
{
public:
    explicit DirichletVectorCorrectionParameters( std::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
    }

    virtual ~DirichletVectorCorrectionParameters() {}

    // This must be a simple variable not a dual or multivariable
    AMP::LinearAlgebra::Variable::shared_ptr d_variable;
};
} // namespace Operator
} // namespace AMP

#endif
