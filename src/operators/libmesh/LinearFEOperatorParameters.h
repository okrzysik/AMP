
#ifndef included_AMP_LinearFEOperatorParameters
#define included_AMP_LinearFEOperatorParameters

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/operators/libmesh/FEOperatorParameters.h"

namespace AMP::Operator {

class LinearFEOperatorParameters : public FEOperatorParameters
{
public:
    /**
      Constructor.
      */
    explicit LinearFEOperatorParameters( std::shared_ptr<AMP::Database> db )
        : FEOperatorParameters( db )
    {
    }

    /**
      Destructor.
      */
    virtual ~LinearFEOperatorParameters() {}

    std::shared_ptr<AMP::Discretization::DOFManager> d_inDofMap;
    std::shared_ptr<AMP::Discretization::DOFManager> d_outDofMap;

protected:
private:
};
} // namespace AMP::Operator

#endif
