#include "CoupledOperator.h"
#include "vectors/Vector.h"
#include "utils/Utilities.h"
#include "utils/ProfilerApp.h"
#include <vector>

namespace AMP {
namespace Operator {


CoupledOperator::CoupledOperator(const boost::shared_ptr<OperatorParameters>& params)
    : ColumnOperator(params)
{
    boost::shared_ptr<CoupledOperatorParameters> myparams = boost::dynamic_pointer_cast<CoupledOperatorParameters>(params);
    d_Operators.push_back(myparams->d_NodeToGaussPointOperator);
    d_Operators.push_back(myparams->d_CopyOperator);
    d_Operators.push_back(myparams->d_MapOperator);
    d_Operators.push_back(myparams->d_BVPOperator);
}


void CoupledOperator::apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
		   const AMP::LinearAlgebra::Vector::shared_ptr &u,
		   AMP::LinearAlgebra::Vector::shared_ptr &r,
		   const double a,
		   const double b)
{
    PROFILE_START("apply");
    if(d_Operators[0]) {
      AMP::LinearAlgebra::Vector::shared_ptr nullVec;
        d_Operators[0]->apply(nullVec,u,d_frozenGaussPointVector,1,0);
    }
    if(d_Operators[1]) {
      if(d_Operators[0]) {
        d_Operators[1]->apply(f,d_frozenGaussPointVector,r,a,b);
      } else {
        d_Operators[1]->apply(f,u,r,a,b);
      }
    }
    if(d_Operators[2]) {
        d_Operators[2]->apply(f,u,r,a,b);
    }
    d_Operators[3]->apply(f,u,r,a,b);
    PROFILE_STOP("apply");
}


}
}


