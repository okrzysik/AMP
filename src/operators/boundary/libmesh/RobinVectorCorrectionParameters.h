
#ifndef included_AMP_RobinVectorCorrectionParameters
#define included_AMP_RobinVectorCorrectionParameters

#include "operators/OperatorParameters.h"

namespace AMP {
namespace Operator {

  class RobinVectorCorrectionParameters : public OperatorParameters {
    public :

      RobinVectorCorrectionParameters(const boost::shared_ptr<AMP::Database> &db)
        : OperatorParameters(db) {  }

      virtual ~RobinVectorCorrectionParameters() { }

      AMP::LinearAlgebra::Variable::shared_ptr d_variable;

  };

}
}

#endif

