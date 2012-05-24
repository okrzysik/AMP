
#ifndef included_AMP_NavierStokesLSWFFEOperatorParameters
#define included_AMP_NavierStokesLSWFFEOperatorParameters

#include "FlowTransportModel.h"
#include "NavierStokesConstants.h"
#include "FEOperatorParameters.h"
#include "vectors/Vector.h"

namespace AMP {
namespace Operator {

  class NavierStokesLSWFFEOperatorParameters : public FEOperatorParameters {
    public :

      NavierStokesLSWFFEOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
        : FEOperatorParameters(db) { }

      virtual ~NavierStokesLSWFFEOperatorParameters() { }

//      boost::shared_ptr<AMP::Discretization::DOFManager> d_dofMap[NavierStokes::TOTAL_NUMBER_OF_VARIABLES];
//      AMP::LinearAlgebra::Vector::shared_ptr d_frozenVec[NavierStokes::TOTAL_NUMBER_OF_VARIABLES]; 
                                              
      boost::shared_ptr<AMP::Discretization::DOFManager> d_dofMap;
      AMP::LinearAlgebra::Vector::shared_ptr d_frozenVec; 

      boost::shared_ptr<FlowTransportModel> d_transportModel;

    protected :

    private :

  };

}
}

#endif

