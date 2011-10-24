
#ifndef included_AMP_FEOperatorParameters
#define included_AMP_FEOperatorParameters

#include "operators/OperatorParameters.h"
#include "ElementOperation.h"

namespace AMP {
namespace Operator {

  /**
   * This class encapsulates parameters used to initialize or reset operators using a
   finite element discretization (FEOperators). It is an abstract base class.    
   @see LinearFEOperator
   @see NonlinearFEOperator
   */
  class FEOperatorParameters : public OperatorParameters {
    public :

      /**
        Constructor.
        */
      FEOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
        : OperatorParameters(db) {  }

      /**
        Destructor.
        */
      virtual ~FEOperatorParameters() {}

      boost::shared_ptr<ElementOperation> d_elemOp; /**< Shared pointer to an element operation */

    protected :

    private :

  };

}
}

#endif

