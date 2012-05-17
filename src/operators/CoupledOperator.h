#ifndef included_AMP_CoupledOperator
#define included_AMP_CoupledOperator

#include "ColumnOperator.h"
#include "CoupledOperatorParameters.h"
#include "vectors/Vector.h"
#include "utils/Utilities.h"
#include <vector>

namespace AMP {
namespace Operator {

class CoupledOperator : public ColumnOperator
{
public :
    CoupledOperator(const boost::shared_ptr<OperatorParameters>& params);

    void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
           const AMP::LinearAlgebra::Vector::shared_ptr &u,
           AMP::LinearAlgebra::Vector::shared_ptr &r,
           const double a = -1.0,
           const double b = 1.0);

      boost::shared_ptr<OperatorParameters>
        getJacobianParameters(const AMP::LinearAlgebra::Vector::shared_ptr & u)
        {
          return (d_Operators[2]->getJacobianParameters(u));
        }

      boost::shared_ptr<AMP::Operator::Operator> getBVPOperator() {
        return d_Operators[2];
      }

      virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
        return d_Operators[2]->getOutputVariable();
      }

      virtual void append(boost::shared_ptr< Operator > op) {
        AMP_ASSERT(d_Operators.size() < 3);
        AMP_ASSERT(op.get() != NULL);
        d_Operators.push_back(op);
      }

      bool isValidInput(boost::shared_ptr<AMP::LinearAlgebra::Vector> &u)
      {
        return d_Operators[2]->isValidInput(u);
      }

      virtual ~CoupledOperator() { }

};


}
}

#endif

