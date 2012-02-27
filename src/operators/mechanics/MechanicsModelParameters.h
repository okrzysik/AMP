
#ifndef included_AMP_MechanicsModelParameters
#define included_AMP_MechanicsModelParameters

#include "operators/ElementPhysicsModelParameters.h"

#include "vectors/Vector.h"

namespace AMP {
  namespace Operator {

    /** A class encapsulating all the parameters that the material model
     * requires to evaluate the stress and/or tangent  */
    class MechanicsModelParameters : public ElementPhysicsModelParameters 
    {
      public :

        /** Constructor */
        MechanicsModelParameters(const boost::shared_ptr<AMP::Database> & db)
          : ElementPhysicsModelParameters(db) { }

        /** Destructor */
        virtual ~MechanicsModelParameters() { }

        /** A vector of deformation gradient values, which are required to
         * compute the stress and/or tangent. */
        boost::shared_ptr<AMP::LinearAlgebra::Vector> d_deformationGradient;

    };

  }
}

#endif


