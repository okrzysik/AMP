
#ifndef included_AMP_NavierStokesGalWFFEOperatorParameters
#define included_AMP_NavierStokesGalWFFEOperatorParameters

#include "FlowTransportModel.h"
#include "FEOperatorParameters.h"
#include "vectors/Vector.h"

namespace AMP {
namespace Operator {

  class NavierStokesGalWFFEOperatorParameters : public FEOperatorParameters {
    public :

      NavierStokesGalWFFEOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
        : FEOperatorParameters(db) { }

      virtual ~NavierStokesGalWFFEOperatorParameters() { }

      AMP::LinearAlgebra::Vector::shared_ptr d_FrozenTemperature; 
                                              
      boost::shared_ptr<FlowTransportModel> d_transportModel;

    protected :

    private :

  };

}
}

#endif

