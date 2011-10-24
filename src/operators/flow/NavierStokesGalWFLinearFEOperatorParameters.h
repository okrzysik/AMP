
#ifndef included_AMP_NavierStokesGalWFLinearFEOperatorParameters
#define included_AMP_NavierStokesGalWFLinearFEOperatorParameters

#include "FEOperatorParameters.h"
#include "FlowTransportModel.h"

#include <vector>

namespace AMP {
namespace Operator {

  class NavierStokesGalWFLinearFEOperatorParameters : public FEOperatorParameters {
    public :

      NavierStokesGalWFLinearFEOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
        : FEOperatorParameters(db) {  }

      virtual ~NavierStokesGalWFLinearFEOperatorParameters() { }

      boost::shared_ptr<FlowTransportModel> d_transportModel; 

      std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_frozenVec; 

    protected :

    private :

  };

}
}

#endif


