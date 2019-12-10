
#ifndef included_AMP_RobinMatrixCorrectionParameters
#define included_AMP_RobinMatrixCorrectionParameters

#include "AMP/operators/boundary/LinearBoundaryOperatorParameters.h"
#include "AMP/operators/boundary/libmesh/NeumannVectorCorrectionParameters.h"
#include "AMP/operators/boundary/libmesh/RobinPhysicsModel.h"

namespace AMP {
namespace Operator {

class RobinMatrixCorrectionParameters : public LinearBoundaryOperatorParameters
{
public:
    explicit RobinMatrixCorrectionParameters( AMP::shared_ptr<AMP::Database> db )
        : LinearBoundaryOperatorParameters( db )
    {
    }

    virtual ~RobinMatrixCorrectionParameters() {}

    AMP::LinearAlgebra::Variable::shared_ptr d_variable;

    AMP::LinearAlgebra::Vector::shared_ptr d_variableFlux;

    std::vector<AMP::LinearAlgebra::Vector::const_shared_ptr> d_elementInputVec;

    AMP::shared_ptr<RobinPhysicsModel> d_robinPhysicsModel;

    AMP::Discretization::DOFManager::shared_ptr d_DofMap;
};
} // namespace Operator
} // namespace AMP

#endif
