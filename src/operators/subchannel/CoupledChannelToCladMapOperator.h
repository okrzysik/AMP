#ifndef included_AMP_CoupledChannelToCladMapOperator
#define included_AMP_CoupledChannelToCladMapOperator

#include "operators/Operator.h"
#include "operators/subchannel/SubchannelPhysicsModel.h"
#include "operators/subchannel/CoupledChannelToCladMapOperatorParameters.h"
#include "vectors/Vector.h"
#include "utils/Utilities.h"
#include <vector>

namespace AMP {
namespace Operator {

  class CoupledChannelToCladMapOperator : public Operator
  {
    public :
      CoupledChannelToCladMapOperator(const boost::shared_ptr<CoupledChannelToCladMapOperatorParameters>& params);

      virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() {
        return d_flowVariable;
      }

      virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
        return d_mapOperator->getOutputVariable();
      }

      void apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u,
          AMP::LinearAlgebra::Vector::shared_ptr r, const double a = -1.0, const double b = 1.0);

      virtual ~CoupledChannelToCladMapOperator() { }
      
    protected :

    private :

      AMP::LinearAlgebra::Variable::shared_ptr d_flowVariable; 

      AMP::LinearAlgebra::Vector::shared_ptr d_subchannelTemperature;

      boost::shared_ptr< AMP::Operator::Operator> d_mapOperator;

      boost::shared_ptr<SubchannelPhysicsModel> d_subchannelPhysicsModel;


  };

}
}

#endif

