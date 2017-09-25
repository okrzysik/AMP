
#ifndef included_AMP_MechanicsNonlinearFEOperatorParameters
#define included_AMP_MechanicsNonlinearFEOperatorParameters

#include "discretization/DOF_Manager.h"
#include "operators/libmesh/FEOperatorParameters.h"
#include "operators/mechanics/MechanicsConstants.h"
#include "operators/mechanics/MechanicsMaterialModel.h"
#include "vectors/Vector.h"

namespace AMP {
namespace Operator {

/**
 * This class encapsulates parameters used to initialize or reset the
 MechanicsNonlinearFEOperator.
 @see MechanicsNonlinearFEOperator
 */
class MechanicsNonlinearFEOperatorParameters : public FEOperatorParameters
{
public:
    /**
      Constructor.
      */
    explicit MechanicsNonlinearFEOperatorParameters( const AMP::shared_ptr<AMP::Database> &db )
        : FEOperatorParameters( db )
    {
    }

    /**
      Destructor.
      */
    virtual ~MechanicsNonlinearFEOperatorParameters() {}

    AMP::shared_ptr<AMP::Discretization::DOFManager> d_dofMap[Mechanics::TOTAL_NUMBER_OF_VARIABLES];

    AMP::shared_ptr<MechanicsMaterialModel> d_materialModel; /**< Material model. */

    AMP::LinearAlgebra::Vector::shared_ptr d_ReferenceTemperature; /**< Reference temperature */

    AMP::LinearAlgebra::Vector::shared_ptr d_EquilibriumVec[Mechanics::TOTAL_NUMBER_OF_VARIABLES];

    AMP::LinearAlgebra::Vector::shared_ptr d_FrozenVec[Mechanics::TOTAL_NUMBER_OF_VARIABLES];

protected:
private:
};
} // namespace Operator
} // namespace AMP

#endif
