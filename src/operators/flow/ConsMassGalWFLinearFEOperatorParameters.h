
#ifndef included_AMP_ConsMassGalWFLinearFEOperatorParameters
#define included_AMP_ConsMassGalWFLinearFEOperatorParameters

#include "LinearFEOperatorParameters.h"
#include "FlowTransportModel.h"

#include <vector>

namespace AMP {
namespace Operator {

  class ConsMassGalWFLinearFEOperatorParameters : public LinearFEOperatorParameters {
    public :

      ConsMassGalWFLinearFEOperatorParameters(const boost::shared_ptr<AMP::Database> &db): LinearFEOperatorParameters(db) {
          d_frozenVec.resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);
      }

      virtual ~ConsMassGalWFLinearFEOperatorParameters() { }

      boost::shared_ptr<FlowTransportModel> d_transportModel; 

      std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_frozenVec; 

    protected :

    private :

  };

}
}

#endif


