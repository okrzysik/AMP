
#ifndef included_AMP_CoupledOperatorParameters
#define included_AMP_CoupledOperatorParameters

/*AMP files */
#include "ColumnOperatorParameters.h"
#include "operators/Operator.h"
/*Boost files */
#include "utils/shared_ptr.h"

#include <vector>

namespace AMP {
namespace Operator {

  /**
    A class that encapsulates the parameters required to construct
    the composite Operator operator.
    @see ColumnOperator
    */
  class CoupledOperatorParameters : public ColumnOperatorParameters {
    public :

      explicit CoupledOperatorParameters(const AMP::shared_ptr<AMP::Database>& db)
        : ColumnOperatorParameters(db) { }

      virtual ~CoupledOperatorParameters() { }

      AMP::shared_ptr<Operator> d_NodeToGaussPointOperator;

      AMP::shared_ptr<Operator> d_CopyOperator;

      AMP::shared_ptr<Operator> d_MapOperator;

      AMP::shared_ptr<Operator> d_BVPOperator;
  };

}  
}


#endif



