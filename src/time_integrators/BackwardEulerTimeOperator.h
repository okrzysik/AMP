#ifndef included_BackwardEulerTimeOperator
#define included_BackwardEulerTimeOperator

#include "TimeOperator.h"

namespace AMP{
namespace TimeIntegrator{

class BackwardEulerTimeOperator: public TimeOperator
{
 public:

  BackwardEulerTimeOperator(AMP::shared_ptr<AMP::Operator::OperatorParameters > params);
  
  void apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, 
         AMP::LinearAlgebra::Vector::const_shared_ptr u,
         AMP::LinearAlgebra::Vector::shared_ptr r,
         const double a = -1.0, const double b=1.0);

  
 protected:
 private:
  
};

}
}

#endif
