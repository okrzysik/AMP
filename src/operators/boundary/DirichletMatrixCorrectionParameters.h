
#ifndef included_AMP_DirichletMatrixCorrectionParameters
#define included_AMP_DirichletMatrixCorrectionParameters

#include "operators/LinearBoundaryOperatorParameters.h"

namespace AMP {
namespace Operator {

  class DirichletMatrixCorrectionParameters : public LinearBoundaryOperatorParameters {
    public :

      DirichletMatrixCorrectionParameters(const boost::shared_ptr<AMP::Database> &db)
        : LinearBoundaryOperatorParameters(db) {  }

      ~DirichletMatrixCorrectionParameters() { }

      //This must be a simple variable not a dual or multivariable
      AMP::LinearAlgebra::Variable::shared_ptr d_variable;

  };

}
}

#endif

