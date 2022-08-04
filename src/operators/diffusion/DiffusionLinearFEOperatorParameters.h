#ifndef included_AMP_DiffusionLinearFEOperatorParameters
#define included_AMP_DiffusionLinearFEOperatorParameters

#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/libmesh/LinearFEOperatorParameters.h"

namespace AMP::Operator {

class DiffusionLinearFEOperatorParameters : public LinearFEOperatorParameters
{
public:
    explicit DiffusionLinearFEOperatorParameters( std::shared_ptr<AMP::Database> db )
        : LinearFEOperatorParameters( db )
    {
    }

    virtual ~DiffusionLinearFEOperatorParameters() {}

    std::shared_ptr<DiffusionTransportModel> d_transportModel;

    std::shared_ptr<AMP::LinearAlgebra::Vector> d_temperature;

    std::shared_ptr<AMP::LinearAlgebra::Vector> d_concentration;

    std::shared_ptr<AMP::LinearAlgebra::Vector> d_burnup;

protected:
private:
};
} // namespace AMP::Operator

#endif
