
#ifndef included_AMP_VolumeIntegralOperatorParameters
#define included_AMP_VolumeIntegralOperatorParameters

#include "FEOperatorParameters.h"
#include "SourcePhysicsModel.h"
#include "discretization/DOF_Manager.h"

namespace AMP {
namespace Operator {

  class VolumeIntegralOperatorParameters : public FEOperatorParameters {
    public :

      VolumeIntegralOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
        : FEOperatorParameters(db) {  }

      ~VolumeIntegralOperatorParameters() { }

      AMP::LinearAlgebra::Vector::shared_ptr d_auxVec;

      boost::shared_ptr<SourcePhysicsModel> d_sourcePhysicsModel;

      AMP::LinearAlgebra::Variable::shared_ptr d_variable;

      AMP::LinearAlgebra::Vector::shared_ptr d_pVector;
    protected :

    private :

  };

}
}

#endif

