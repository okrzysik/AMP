#ifndef included_AMP_VectorCopyOperator
#define  included_AMP_VectorCopyOperator

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

#include "boost/shared_ptr.hpp"
#include "operators/VectorCopyOperatorParameters.h"
#include "operators/Operator.h"
#include "vectors/Vector.h"

namespace AMP {
namespace Operator {

class VectorCopyOperator : public Operator {
public:

  VectorCopyOperator(const boost::shared_ptr<VectorCopyOperatorParameters> &params);

  virtual ~VectorCopyOperator(){}
  
  virtual void apply(AMP::LinearAlgebra::Vector::const_shared_ptr f,
		     AMP::LinearAlgebra::Vector::const_shared_ptr u,
		     AMP::LinearAlgebra::Vector::shared_ptr r,
		     const double a = -1.0,
		     const double b = 1.0);
  
  AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable();

  AMP::LinearAlgebra::Variable::shared_ptr getInputVariable();

 private:
  // vector to copy into
  boost::shared_ptr<AMP::LinearAlgebra::Vector> d_copyVector;
  boost::shared_ptr<AMP::LinearAlgebra::Variable> d_copyVariable;
  
};
 
} // namespace Operator
} // namespace AMP

#endif
