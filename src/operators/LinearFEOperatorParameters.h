
#ifndef included_AMP_LinearFEOperatorParameters
#define included_AMP_LinearFEOperatorParameters

#include "operators/FEOperatorParameters.h"
#include "discretization/DOF_Manager.h"

namespace AMP {
  namespace Operator {

   class LinearFEOperatorParameters : public FEOperatorParameters {
      public :

        /**
          Constructor.
          */
        LinearFEOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
          : FEOperatorParameters(db) {  }

        /**
          Destructor.
          */
        virtual ~LinearFEOperatorParameters() {}

        boost::shared_ptr<AMP::Discretization::DOFManager> d_inDofMap;
        boost::shared_ptr<AMP::Discretization::DOFManager> d_outDofMap;

      protected :

      private :

    };

  }
}

#endif

