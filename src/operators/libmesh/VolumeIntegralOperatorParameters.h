
#ifndef included_AMP_VolumeIntegralOperatorParameters
#define included_AMP_VolumeIntegralOperatorParameters

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/operators/libmesh/FEOperatorParameters.h"
#include "AMP/operators/libmesh/SourcePhysicsModel.h"

namespace AMP::Operator {


class VolumeIntegralOperatorParameters : public FEOperatorParameters
{
public:
    explicit VolumeIntegralOperatorParameters( std::shared_ptr<AMP::Database> db )
        : FEOperatorParameters( db )
    {
    }

    virtual ~VolumeIntegralOperatorParameters() {}

    AMP::LinearAlgebra::Vector::shared_ptr d_auxVec;
    std::shared_ptr<SourcePhysicsModel> d_sourcePhysicsModel;
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_variable;
    AMP::LinearAlgebra::Vector::shared_ptr d_pVector;

protected:
private:
};


} // namespace AMP::Operator

#endif
