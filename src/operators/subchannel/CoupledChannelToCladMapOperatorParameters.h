
#ifndef included_AMP_CoupledChannelToCladMapOperatorParameters
#define included_AMP_CoupledChannelToCladMapOperatorParameters

/*AMP files */
#include "operators/OperatorParameters.h"
#include "operators/subchannel/SubchannelPhysicsModel.h"
#include "operators/Operator.h"
/*Boost files */
#include "boost/shared_ptr.hpp"

#include <vector>

namespace AMP {
namespace Operator {

  /**
    A class that encapsulates the parameters required to construct
    the composite Operator operator.
    @see ColumnOperator
    */
  class CoupledChannelToCladMapOperatorParameters : public OperatorParameters {
    public :

      CoupledChannelToCladMapOperatorParameters(const boost::shared_ptr<AMP::Database>& db)
        : OperatorParameters(db) { }

      virtual ~CoupledChannelToCladMapOperatorParameters() { }

      boost::shared_ptr<AMP::Operator::Operator> d_mapOperator ; 
      boost::shared_ptr<AMP::LinearAlgebra::Variable> d_variable ; 
      boost::shared_ptr<AMP::Mesh::Mesh> d_subchannelMesh ; 
      boost::shared_ptr<AMP::LinearAlgebra::Vector> d_vector ; 
      boost::shared_ptr<AMP::Operator::SubchannelPhysicsModel> d_subchannelPhysicsModel; 

  };

}  
}


#endif



