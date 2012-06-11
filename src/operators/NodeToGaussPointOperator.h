#ifndef included_AMP_NodeToGaussPointOperator
#define included_AMP_NodeToGaussPointOperator

#include "operators/Operator.h"

namespace AMP {
namespace Operator {


class NodeToGaussPointOperator : public Operator
{
public :

  NodeToGaussPointOperator (const boost::shared_ptr<OperatorParameters> & params) : Operator (params) {
    d_NodalVariable.reset(new AMP::LinearAlgebra::Variable(params->d_db->getString("InputVariable")));
    d_GaussPtVariable.reset(new AMP::LinearAlgebra::Variable(params->d_db->getString("OutputVariable")));
  }

  virtual ~NodeToGaussPointOperator() { }

  void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
        AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);
    
protected :
  AMP::LinearAlgebra::Variable::shared_ptr d_NodalVariable;
  AMP::LinearAlgebra::Variable::shared_ptr d_GaussPtVariable;


};

}
}

#endif
