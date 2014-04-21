
#ifndef included_AMP_NavierStokesLinearFEOperatorParameters
#define included_AMP_NavierStokesLinearFEOperatorParameters

#include "operators/flow/FlowTransportModel.h"
#include "operators/flow/NavierStokesConstants.h"
#include "operators/libmesh/LinearFEOperatorParameters.h"
#include "vectors/Vector.h"

namespace AMP {
namespace Operator {

  class NavierStokesLinearFEOperatorParameters : public LinearFEOperatorParameters {
    public :

      NavierStokesLinearFEOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
        : LinearFEOperatorParameters(db) { }

      virtual ~NavierStokesLinearFEOperatorParameters() { }

//      AMP::LinearAlgebra::Vector::shared_ptr d_frozenVec[NavierStokes::TOTAL_NUMBER_OF_VARIABLES]; 
      AMP::LinearAlgebra::Vector::shared_ptr d_frozenVec; 
                                              
      boost::shared_ptr<FlowTransportModel> d_transportModel;

    protected :

    private :

  };

}
}

#endif

