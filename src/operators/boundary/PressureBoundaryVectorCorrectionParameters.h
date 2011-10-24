
#ifndef included_AMP_PressureBoundaryVectorCorrectionParameters
#define included_AMP_PressureBoundaryVectorCorrectionParameters

#include "operators/LinearBoundaryOperatorParameters.h"
#include "operators/OperatorParameters.h"

namespace AMP {
namespace Operator {

  class PressureBoundaryVectorCorrectionParameters : public OperatorParameters {
    public :

      PressureBoundaryVectorCorrectionParameters(const boost::shared_ptr<AMP::Database> &db)
        : OperatorParameters(db) {  }

      ~PressureBoundaryVectorCorrectionParameters() { }

      AMP::LinearAlgebra::Variable::shared_ptr d_variable;
 
      AMP::LinearAlgebra::Vector::shared_ptr d_variablePressure;

  };

}
}

#endif

