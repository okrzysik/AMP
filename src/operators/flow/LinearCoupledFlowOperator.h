#ifndef included_AMP_LinearCoupledFlowOperator
#define included_AMP_LinearCoupledFlowOperator

#include "LinearCoupledFlowOperatorParameters.h"
#include "vectors/Vector.h"
#include "utils/Utilities.h"
#include <vector>

namespace AMP {
namespace Operator {

  class LinearCoupledFlowOperator : public Operator
  {
    public :
      LinearCoupledFlowOperator(const boost::shared_ptr<OperatorParameters>& params)
        : Operator(params){ (void) params;
        }

      virtual ~LinearCoupledFlowOperator() { }

      virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u, 
          AMP::LinearAlgebra::Vector::shared_ptr &r, const double a = -1.0, const double b = 1.0);

      virtual void reset(const boost::shared_ptr<OperatorParameters>& params);

      virtual void append(boost::shared_ptr< Operator > op);

      virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() ;

      virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1);

    protected :

      std::vector< boost::shared_ptr< Operator > > d_Operators;

    private :

  };

}
}

#endif

