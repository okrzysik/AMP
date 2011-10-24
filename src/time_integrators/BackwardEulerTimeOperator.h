#ifndef included_BackwardEulerTimeOperator
#define included_BackwardEulerTimeOperator

#include "TimeOperator.h"

namespace AMP{
namespace TimeIntegrator{

class BackwardEulerTimeOperator: public TimeOperator
{
 public:

  BackwardEulerTimeOperator(boost::shared_ptr<AMP::Operator::OperatorParameters > params);
  
  void apply(const boost::shared_ptr<AMP::LinearAlgebra::Vector>  &f, 
         const boost::shared_ptr<AMP::LinearAlgebra::Vector>  &u,
         boost::shared_ptr<AMP::LinearAlgebra::Vector>  &r,
         const double a = -1.0, const double b=1.0);

  
 protected:
 private:
  
};

}
}

#endif
