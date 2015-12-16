
#ifndef included_AMP_NavierStokesGalWFFEOperatorParameters
#define included_AMP_NavierStokesGalWFFEOperatorParameters

#include "operators/flow/FlowTransportModel.h"
#include "operators/libmesh/LinearFEOperatorParameters.h"
#include "vectors/Vector.h"

namespace AMP {
namespace Operator {

class NavierStokesGalWFFEOperatorParameters : public LinearFEOperatorParameters {
public:
    explicit NavierStokesGalWFFEOperatorParameters( const AMP::shared_ptr<AMP::Database> &db )
        : LinearFEOperatorParameters( db )
    {
    }

    virtual ~NavierStokesGalWFFEOperatorParameters() {}

    AMP::LinearAlgebra::Vector::shared_ptr d_FrozenTemperature;

    AMP::shared_ptr<FlowTransportModel> d_transportModel;

protected:
private:
};
}
}

#endif
