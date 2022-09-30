#ifndef included_AMP_DiffusionNonlinearFEOperatorParameters
#define included_AMP_DiffusionNonlinearFEOperatorParameters

#include <vector>

#include "AMP/operators/diffusion/DiffusionConstants.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/libmesh/FEOperatorParameters.h"
#include "AMP/vectors/Vector.h"


namespace AMP::Operator {

class DiffusionNonlinearFEOperatorParameters : public FEOperatorParameters
{
public:
    explicit DiffusionNonlinearFEOperatorParameters( std::shared_ptr<AMP::Database> db )
        : FEOperatorParameters( db )
    {
    }

    virtual ~DiffusionNonlinearFEOperatorParameters() {}

    std::shared_ptr<DiffusionTransportModel> d_transportModel;

    std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> d_FrozenVecs;
};

} // namespace AMP::Operator

#endif
