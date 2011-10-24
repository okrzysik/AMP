
#ifndef included_AMP_ConsMomentumGalWFFEOperatorParameters
#define included_AMP_ConsMomentumGalWFFEOperatorParameters

#include "FlowTransportModel.h"
#include "FEOperatorParameters.h"
#include "vectors/Vector.h"

namespace AMP {
namespace Operator {

  class ConsMomentumGalWFFEOperatorParameters : public FEOperatorParameters {
    public :

      ConsMomentumGalWFFEOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
        : FEOperatorParameters(db) { }

      virtual ~ConsMomentumGalWFFEOperatorParameters() { }

      AMP::LinearAlgebra::Vector::shared_ptr d_FrozenTemperature; 
                                              
      AMP::LinearAlgebra::Vector::shared_ptr d_FrozenPressure; 
                                              
      boost::shared_ptr<FlowTransportModel> d_transportModel;

    protected :

    private :

  };

}
}

#endif

