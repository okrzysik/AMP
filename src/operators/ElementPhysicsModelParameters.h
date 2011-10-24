
#ifndef included_AMP_ElementPhysicsModelParameters
#define included_AMP_ElementPhysicsModelParameters

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"

#include "utils/ParameterBase.h"

namespace AMP {
namespace Operator {

  /**
    An abstract base class that encapsulates the parameters used to
    initialize the ElementPhysicsModel (Constitutive/Material model)
    used within a FEOperator. 
    @see ElementPhysicsModel
    */
  class ElementPhysicsModelParameters : public ParameterBase
  {
    public :

      /** 
        Constructor.
        */
      ElementPhysicsModelParameters(const boost::shared_ptr<AMP::Database> & db)
        : d_db(db) {}

      /**
        Destructor.
        */
      virtual ~ElementPhysicsModelParameters() { }

      boost::shared_ptr<AMP::Database> d_db; /**< Database object which needs to be 
                                               initialized specific to the material model.  */

    protected :

    private :

  };

}
}

#endif


