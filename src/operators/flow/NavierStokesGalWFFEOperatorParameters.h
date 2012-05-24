
#ifndef included_AMP_NavierStokesGalWFFEOperatorParameters
#define included_AMP_NavierStokesGalWFFEOperatorParameters

#include "FlowTransportModel.h"
#include "LinearFEOperatorParameters.h"
#include "vectors/Vector.h"

namespace AMP {
namespace Operator {

  class NavierStokesGalWFFEOperatorParameters : public LinearFEOperatorParameters {
    public :

      NavierStokesGalWFFEOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
        : LinearFEOperatorParameters(db) { }

      virtual ~NavierStokesGalWFFEOperatorParameters() { }

      AMP::LinearAlgebra::Vector::shared_ptr d_FrozenTemperature; 
                                              
      boost::shared_ptr<FlowTransportModel> d_transportModel;

    protected :

    private :

  };

}
}

#endif

