#include "operators/VectorCopyOperator.h"

namespace AMP{
namespace Operator{

VectorCopyOperator::VectorCopyOperator(const boost::shared_ptr<VectorCopyOperatorParameters> &params):AMP::Operator::Operator(params)
{
  boost::shared_ptr<AMP::Operator::VectorCopyOperatorParameters> copyParams = boost::dynamic_pointer_cast< AMP::Operator::VectorCopyOperatorParameters>(params);
  d_copyVariable = copyParams->d_copyVariable;
  d_copyVector = copyParams->d_copyVector;
  
  AMP_INSIST(d_copyVector, "must have non NULL CopyVector");
  AMP_INSIST(d_copyVariable, "must have non NULL CopyVeriable");
   
}

void
VectorCopyOperator::apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
			  const  AMP::LinearAlgebra::Vector::shared_ptr &u,
			  AMP::LinearAlgebra::Vector::shared_ptr  &r,
			  const double a,
			  const double b)
{
  AMP::LinearAlgebra::Vector::shared_ptr vecToCopy = u->subsetVectorForVariable( this->getOutputVariable() );
  d_copyVector->copyVector(vecToCopy);
  
}

AMP::LinearAlgebra::Variable::shared_ptr
VectorCopyOperator::getOutputVariable()
{
  return d_copyVariable;
}

AMP::LinearAlgebra::Variable::shared_ptr
VectorCopyOperator::getInputVariable()
{
  return d_copyVariable;
}

  
}
}
