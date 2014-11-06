#include "BackwardEulerTimeOperator.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

namespace AMP{
namespace TimeIntegrator{

BackwardEulerTimeOperator::BackwardEulerTimeOperator(AMP::shared_ptr<AMP::Operator::OperatorParameters > params):TimeOperator(params)
{
}

void
BackwardEulerTimeOperator::apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, 
                 AMP::LinearAlgebra::Vector::const_shared_ptr u,
                 AMP::LinearAlgebra::Vector::shared_ptr r,
                 const double a, const double b)
{

  // this routine evaluates a*[ ( M(u)-M(uOld) )/dt-fRhs(u) -source_term] +b*f
  // where the time operator is given by u_t = fRhs(u) 

  AMP::shared_ptr<AMP::LinearAlgebra::Vector>  fTmp;

  AMP_INSIST(d_pRhsOperator.get()!=NULL, "ERROR: AMP::TimeIntegrator::TimeIntegrator::TimeOperator::apply, the rhs operator is NULL!");
  
  if(d_pMassOperator.get()!=NULL)
    {
      if(d_bLinearMassOperator)
    {
      d_pScratchVector->subtract(*u, *d_pPreviousTimeSolution);
      d_pMassOperator->apply(fTmp,d_pScratchVector,r,1.0/d_dCurrentDt,0.0);
    }
      else
    {
      d_pMassOperator->apply(fTmp, d_pPreviousTimeSolution,d_pScratchVector,1.0,0.0);
      d_pMassOperator->apply(fTmp, u, r,-1.0,0.0);
      r->add(*r, *d_pScratchVector);
      r->scale(1.0/d_dCurrentDt);
    }
    }
  else
    {
      r->subtract(*u, *d_pPreviousTimeSolution);
      r->scale(1.0/d_dCurrentDt);      
    }
  
  d_pRhsOperator->apply(fTmp, u, d_pScratchVector);

  r->add(*r, *d_pScratchVector);

  //subtract any source sink terms and boundary corrections

  r->subtract(*r, *d_pSourceTerm);
  
  if(f.get()==NULL)
    {
      r->scale(a);
    }
  else
    {
      r->axpby(b, a, *f);
    }

}
  
}
}

