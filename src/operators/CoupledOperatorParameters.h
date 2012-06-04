
#ifndef included_AMP_CoupledOperatorParameters
#define included_AMP_CoupledOperatorParameters

/*AMP files */
#include "ColumnOperatorParameters.h"
#include "operators/Operator.h"
/*Boost files */
#include "boost/shared_ptr.hpp"

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

      CoupledOperatorParameters(const boost::shared_ptr<AMP::Database>& db)
        : ColumnOperatorParameters(db) { }

      virtual ~CoupledOperatorParameters() { }

      boost::shared_ptr<Operator> d_NodeToGaussPointOperator;

      boost::shared_ptr<Operator> d_CopyOperator;

      boost::shared_ptr<Operator> d_MapOperator;

      boost::shared_ptr<Operator> d_BVPOperator;
  };

}  
}


#endif



