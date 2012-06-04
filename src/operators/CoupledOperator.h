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
          return (d_Operators[3]->getJacobianParameters(u));
        }

      boost::shared_ptr<AMP::Operator::Operator> getBVPOperator() {
        return d_Operators[3];
      }

      virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
        return d_Operators[3]->getOutputVariable();
      }

      virtual void append(boost::shared_ptr< Operator > op) {
        AMP_ASSERT(d_Operators.size() < 4);
        AMP_ASSERT(op.get() != NULL);
        d_Operators.push_back(op);
      }

      bool isValidInput(boost::shared_ptr<AMP::LinearAlgebra::Vector> &u)
      {
        return d_Operators[3]->isValidInput(u);
      }

      void setFrozenGaussPointVector(AMP::LinearAlgebra::Vector::shared_ptr u) {
        d_frozenGaussPointVector = u;
      }

      virtual ~CoupledOperator() { }

protected:
      AMP::LinearAlgebra::Vector::shared_ptr d_frozenGaussPointVector;

};


}
}

#endif

