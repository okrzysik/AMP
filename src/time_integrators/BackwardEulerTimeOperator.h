#ifndef included_BackwardEulerTimeOperator
#define included_BackwardEulerTimeOperator

#include "TimeOperator.h"

namespace AMP{
namespace TimeIntegrator{

class BackwardEulerTimeOperator: public TimeOperator
{
 public:

  explicit BackwardEulerTimeOperator(AMP::shared_ptr<AMP::Operator::OperatorParameters > params);
  
  void apply(AMP::LinearAlgebra::Vector::const_shared_ptr u,
	     AMP::LinearAlgebra::Vector::shared_ptr f ) override;

  
 protected:
 private:
  
};

}
}

#endif
