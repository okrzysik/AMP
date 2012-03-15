
#ifndef included_AMP_MechanicsLinearFEOperatorParameters
#define included_AMP_MechanicsLinearFEOperatorParameters

#include "operators/LinearFEOperatorParameters.h"
#include "operators/mechanics/MechanicsMaterialModel.h"

#include <vector>

namespace AMP {
  namespace Operator {

    /**
     * This class encapsulates parameters used to initialize or reset the
     MechanicsLinearFEOperator.    
     @see MechanicsLinearFEOperator
     */
    class MechanicsLinearFEOperatorParameters : public LinearFEOperatorParameters {
      public :

        /**
          Constructor.
          */
        MechanicsLinearFEOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
          : LinearFEOperatorParameters(db) {  }

        /**
          Destructor.
          */
        virtual ~MechanicsLinearFEOperatorParameters() { }

        boost::shared_ptr<MechanicsMaterialModel> d_materialModel; /**< Material model. */

        AMP::LinearAlgebra::Vector::shared_ptr d_dispVec; /**< Displacement vector, which is passed from 
                                                            MechanicsNonlinearFEOperator to MechanicsLinearFEOperator. */

      protected :

      private :

    };

  }
}

#endif


