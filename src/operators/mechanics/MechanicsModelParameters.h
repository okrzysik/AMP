
#ifndef included_AMP_MechanicsModelParameters
#define included_AMP_MechanicsModelParameters

#include "operators/ElementPhysicsModelParameters.h"

#include "vectors/Vector.h"

namespace AMP {
  namespace Operator {

    class MechanicsModelParameters : public ElementPhysicsModelParameters 
    {
      public :

        MechanicsModelParameters(const boost::shared_ptr<AMP::Database> & db)
          : ElementPhysicsModelParameters(db) { }

        virtual ~MechanicsModelParameters() { }

        boost::shared_ptr<AMP::LinearAlgebra::Vector> d_deformationGradient;

    };

  }
}

#endif


