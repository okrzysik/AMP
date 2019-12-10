#ifndef included_AMP_DiffusionLinearFEOperatorParameters
#define included_AMP_DiffusionLinearFEOperatorParameters

#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/libmesh/LinearFEOperatorParameters.h"

namespace AMP {
namespace Operator {

class DiffusionLinearFEOperatorParameters : public LinearFEOperatorParameters
{
public:
    explicit DiffusionLinearFEOperatorParameters( AMP::shared_ptr<AMP::Database> db )
        : LinearFEOperatorParameters( db )
    {
    }

    virtual ~DiffusionLinearFEOperatorParameters() {}

    AMP::shared_ptr<DiffusionTransportModel> d_transportModel;

    AMP::LinearAlgebra::Vector::shared_ptr d_temperature;

    AMP::LinearAlgebra::Vector::shared_ptr d_concentration;

    AMP::LinearAlgebra::Vector::shared_ptr d_burnup;

protected:
private:
};
} // namespace Operator
} // namespace AMP

#endif
