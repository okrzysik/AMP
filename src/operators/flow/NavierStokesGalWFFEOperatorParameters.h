
#ifndef included_AMP_NavierStokesGalWFFEOperatorParameters
#define included_AMP_NavierStokesGalWFFEOperatorParameters

#include "AMP/operators/flow/FlowTransportModel.h"
#include "AMP/operators/libmesh/LinearFEOperatorParameters.h"
#include "AMP/vectors/Vector.h"

namespace AMP {
namespace Operator {

class NavierStokesGalWFFEOperatorParameters : public LinearFEOperatorParameters
{
public:
    explicit NavierStokesGalWFFEOperatorParameters( std::shared_ptr<AMP::Database> db )
        : LinearFEOperatorParameters( db )
    {
    }

    virtual ~NavierStokesGalWFFEOperatorParameters() {}

    AMP::LinearAlgebra::Vector::shared_ptr d_FrozenTemperature;

    std::shared_ptr<FlowTransportModel> d_transportModel;

protected:
private:
};
} // namespace Operator
} // namespace AMP

#endif
