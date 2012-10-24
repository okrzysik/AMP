
#ifndef included_AMP_FlowFrapconJacobianParameters
#define included_AMP_FlowFrapconJacobianParameters

#include "operators/OperatorParameters.h"

namespace AMP {
namespace Operator {

  class FlowFrapconJacobianParameters : public OperatorParameters {
    public :

      FlowFrapconJacobianParameters(const boost::shared_ptr<AMP::Database> &db)
        : OperatorParameters(db) {  }

      virtual ~FlowFrapconJacobianParameters() { }

      AMP::LinearAlgebra::Variable::shared_ptr d_variable;

      AMP::LinearAlgebra::Vector::shared_ptr d_frozenSolution;

  };

}
}

#endif

