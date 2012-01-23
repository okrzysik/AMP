
#ifndef included_AMP_MechanicsNonlinearFEOperatorParameters
#define included_AMP_MechanicsNonlinearFEOperatorParameters

#include "FEOperatorParameters.h"
#include "MechanicsMaterialModel.h"
#include "discretization/DOF_Manager.h"
#include "MechanicsConstants.h"

#include "vectors/Vector.h"

namespace AMP {
  namespace Operator {

    /**
     * This class encapsulates parameters used to initialize or reset the
     MechanicsNonlinearFEOperator.    
     @see MechanicsNonlinearFEOperator
     */
    class MechanicsNonlinearFEOperatorParameters : public FEOperatorParameters {
      public :

        /**
          Constructor.
          */
        MechanicsNonlinearFEOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
          : FEOperatorParameters(db) { }

        /**
          Destructor.
          */
        virtual ~MechanicsNonlinearFEOperatorParameters() { }

        boost::shared_ptr<AMP::Discretization::DOFManager> d_dofMap[Mechanics::TOTAL_NUMBER_OF_VARIABLES];

        boost::shared_ptr<MechanicsMaterialModel> d_materialModel; /**< Material model. */

        AMP::LinearAlgebra::Vector::shared_ptr d_EquilibriumDisplacement; /**< Displacement at equilibrium. */

        AMP::LinearAlgebra::Vector::shared_ptr d_EquilibriumTemperature; /**< Temperature at equilibrium. */

        AMP::LinearAlgebra::Vector::shared_ptr d_EquilibriumBurnup; /**< Burnup at equilibrium. */

        AMP::LinearAlgebra::Vector::shared_ptr d_EquilibriumOxygenConcentration; /**< Oxygen concentration at equilibrium. */

        AMP::LinearAlgebra::Vector::shared_ptr d_EquilibriumLHGR; /**< Linear Heat Generation Rate at equilibrium. */

        AMP::LinearAlgebra::Vector::shared_ptr d_ReferenceTemperature; /**< Reference temperature */

        AMP::LinearAlgebra::Vector::shared_ptr d_FrozenTemperature; /**< Frozen temperature. This is used when displacement 
                                                                      and temperature are solved in an uncoupled fashion. */

        AMP::LinearAlgebra::Vector::shared_ptr d_FrozenBurnup; /**< Frozen burnup. This is used when displacement and burnup
                                                                 are solved in an uncoupled fashion. */

        AMP::LinearAlgebra::Vector::shared_ptr d_FrozenOxygenConcentration; /**< Frozen oxygen concentration. This is used when
                                                                              displacement and oxygen concentration
                                                                              are solved in an uncoupled fashion. */

        AMP::LinearAlgebra::Vector::shared_ptr d_FrozenLHGR; /**< Frozen Linear Heat Generation Rate. This is used when creep is simulated. */

      protected :

      private :

    };

  }
}

#endif

