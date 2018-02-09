
#ifndef included_AMP_NavierStokesLSWFFEOperatorParameters
#define included_AMP_NavierStokesLSWFFEOperatorParameters

#include "AMP/operators/flow/FlowTransportModel.h"
#include "AMP/operators/flow/NavierStokesConstants.h"
#include "AMP/operators/libmesh/FEOperatorParameters.h"
#include "AMP/vectors/Vector.h"

namespace AMP {
namespace Operator {

class NavierStokesLSWFFEOperatorParameters : public FEOperatorParameters
{
public:
    explicit NavierStokesLSWFFEOperatorParameters( const AMP::shared_ptr<AMP::Database> &db )
        : FEOperatorParameters( db )
    {
    }

    virtual ~NavierStokesLSWFFEOperatorParameters() {}

    //      AMP::shared_ptr<AMP::Discretization::DOFManager>
    //      d_dofMap[NavierStokes::TOTAL_NUMBER_OF_VARIABLES];
    //      AMP::LinearAlgebra::Vector::shared_ptr
    //      d_frozenVec[NavierStokes::TOTAL_NUMBER_OF_VARIABLES];

    AMP::shared_ptr<AMP::Discretization::DOFManager> d_dofMap;
    AMP::LinearAlgebra::Vector::shared_ptr d_frozenVec;

    AMP::shared_ptr<FlowTransportModel> d_transportModel;

protected:
private:
};
} // namespace Operator
} // namespace AMP

#endif
