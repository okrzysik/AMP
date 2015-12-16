
#ifndef included_AMP_VolumeIntegralOperatorParameters
#define included_AMP_VolumeIntegralOperatorParameters

#include "discretization/DOF_Manager.h"
#include "operators/libmesh/FEOperatorParameters.h"
#include "operators/libmesh/SourcePhysicsModel.h"

namespace AMP {
namespace Operator {


class VolumeIntegralOperatorParameters : public FEOperatorParameters
{
public:
    explicit VolumeIntegralOperatorParameters( const AMP::shared_ptr<AMP::Database> &db )
        : FEOperatorParameters( db )
    {
    }

    virtual ~VolumeIntegralOperatorParameters() {}

    AMP::LinearAlgebra::Vector::shared_ptr d_auxVec;
    AMP::shared_ptr<SourcePhysicsModel> d_sourcePhysicsModel;
    AMP::LinearAlgebra::Variable::shared_ptr d_variable;
    AMP::LinearAlgebra::Vector::shared_ptr d_pVector;

protected:
private:
};


} // namespace Operator
} // namespace AMP

#endif
