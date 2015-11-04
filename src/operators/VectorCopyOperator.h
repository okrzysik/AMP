#ifndef included_AMP_VectorCopyOperator
#define  included_AMP_VectorCopyOperator

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

#include "utils/shared_ptr.h"
#include "operators/VectorCopyOperatorParameters.h"
#include "operators/Operator.h"
#include "vectors/Vector.h"

namespace AMP {
namespace Operator {

class VectorCopyOperator : public Operator {
public:

  VectorCopyOperator(const AMP::shared_ptr<VectorCopyOperatorParameters> &params);

  virtual ~VectorCopyOperator(){}
  
  virtual void apply(AMP::LinearAlgebra::Vector::const_shared_ptr u,
		     AMP::LinearAlgebra::Vector::shared_ptr f ) override;
  
  AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable();

  AMP::LinearAlgebra::Variable::shared_ptr getInputVariable();

 private:
  // vector to copy into
  AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_copyVector;
  AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_copyVariable;
  
};
 
} // namespace Operator
} // namespace AMP

#endif
