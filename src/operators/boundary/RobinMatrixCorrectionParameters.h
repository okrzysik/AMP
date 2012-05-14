
#ifndef included_AMP_RobinMatrixCorrectionParameters
#define included_AMP_RobinMatrixCorrectionParameters

#include "NeumannVectorCorrectionParameters.h"
#include "operators/boundary/LinearBoundaryOperatorParameters.h"
#include "RobinPhysicsModel.h"

namespace AMP {
namespace Operator {

  class RobinMatrixCorrectionParameters : public LinearBoundaryOperatorParameters {
    public :

      RobinMatrixCorrectionParameters(const boost::shared_ptr<AMP::Database> &db)
        : LinearBoundaryOperatorParameters(db) {  }

      ~RobinMatrixCorrectionParameters() { }

      AMP::LinearAlgebra::Variable::shared_ptr d_variable;

      AMP::LinearAlgebra::Vector::shared_ptr d_variableFlux;
      
      std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_elementInputVec;
      
      boost::shared_ptr<RobinPhysicsModel> d_robinPhysicsModel;

      AMP::Discretization::DOFManager::shared_ptr d_DofMap; 
  };

}
}

#endif

