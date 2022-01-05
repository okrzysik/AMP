#ifndef included_AMP_MassLinearFEOperatorParameters
#define included_AMP_MassLinearFEOperatorParameters

#include "AMP/operators/libmesh/LinearFEOperatorParameters.h"

#include "AMP/operators/libmesh/MassDensityModel.h"

namespace AMP::Operator {

class MassLinearFEOperatorParameters : public LinearFEOperatorParameters
{
public:
    explicit MassLinearFEOperatorParameters( std::shared_ptr<AMP::Database> db )
        : LinearFEOperatorParameters( db )
    {
    }

    virtual ~MassLinearFEOperatorParameters() {}

    std::shared_ptr<MassDensityModel> d_densityModel;

protected:
private:
};
} // namespace AMP::Operator

#endif
