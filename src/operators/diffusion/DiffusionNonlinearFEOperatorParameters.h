#ifndef included_AMP_DiffusionNonlinearFEOperatorParameters
#define included_AMP_DiffusionNonlinearFEOperatorParameters

#include <vector>

#include "operators/diffusion/DiffusionConstants.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/libmesh/FEOperatorParameters.h"
#include "vectors/Vector.h"


namespace AMP {
namespace Operator {

class DiffusionNonlinearFEOperatorParameters : public FEOperatorParameters
{
public:
    explicit DiffusionNonlinearFEOperatorParameters( const AMP::shared_ptr<AMP::Database> &db )
        : FEOperatorParameters( db )
    {
    }

    virtual ~DiffusionNonlinearFEOperatorParameters() {}

    AMP::shared_ptr<DiffusionTransportModel> d_transportModel;

    AMP::LinearAlgebra::Vector::shared_ptr d_FrozenTemperature;

    AMP::LinearAlgebra::Vector::shared_ptr d_FrozenConcentration;

    AMP::LinearAlgebra::Vector::shared_ptr d_FrozenBurnup;

protected:
private:
};
} // namespace Operator
} // namespace AMP

#endif
