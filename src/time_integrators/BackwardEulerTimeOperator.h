#ifndef included_BackwardEulerTimeOperator
#define included_BackwardEulerTimeOperator

#include "TimeOperator.h"

namespace AMP{
namespace TimeIntegrator{

class BackwardEulerTimeOperator: public TimeOperator
{
 public:

  BackwardEulerTimeOperator(boost::shared_ptr<AMP::Operator::OperatorParameters > params);
  
  void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, 
         const AMP::LinearAlgebra::Vector::shared_ptr &u,
         AMP::LinearAlgebra::Vector::shared_ptr &r,
         const double a = -1.0, const double b=1.0);

  
 protected:
 private:
  
};

}
}

#endif
