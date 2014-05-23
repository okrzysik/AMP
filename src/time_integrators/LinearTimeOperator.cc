#include "LinearTimeOperator.h"
#include "TimeOperatorParameters.h"
#include "utils/InputDatabase.h"



namespace AMP{
namespace TimeIntegrator{

LinearTimeOperator::LinearTimeOperator(boost::shared_ptr<AMP::Operator::OperatorParameters > in_params):LinearOperator(in_params)
{
  d_bModifyRhsOperatorMatrix=false;
  
  boost::shared_ptr<TimeOperatorParameters> params = boost::dynamic_pointer_cast<TimeOperatorParameters>(in_params);

  d_pRhsOperator = boost::dynamic_pointer_cast<LinearOperator>(params->d_pRhsOperator);
  d_pMassOperator = boost::dynamic_pointer_cast<LinearOperator>(params->d_pMassOperator);


  d_beta = 1.0;
  reset(in_params);

}

LinearTimeOperator::~LinearTimeOperator()
{
}

void
LinearTimeOperator::getFromInput(const boost::shared_ptr<AMP::Database> &db)
{
    
  AMP_INSIST(db->keyExists("ScalingFactor"), "key ScalingFactor missing in input");

  d_dScalingFactor = db->getDouble("ScalingFactor");
  d_current_time = db->getDoubleWithDefault("CurrentTime",0.0);
  d_bAlgebraicComponent = db->getBoolWithDefault("bAlgebraicComponent", false);

}
  
void
LinearTimeOperator::reset(const boost::shared_ptr<AMP::Operator::OperatorParameters>& in_params)
{
  boost::shared_ptr<TimeOperatorParameters> params = boost::dynamic_pointer_cast<TimeOperatorParameters>(in_params);

  getFromInput(params->d_db);

  // the parameter object for the rhs operator will be NULL during construction, but should
  // not be NULL during a reset based on getJacobianParameters
  if(params->d_pRhsOperatorParameters.get()!=NULL)
    {
      d_pRhsOperator->reset(params->d_pRhsOperatorParameters);
    }

  // the parameter object for the mass operator will be NULL during construction, but should
  // not be NULL during a reset based on getJacobianParameters
  if(params->d_pMassOperatorParameters.get()!=NULL)
    {
      d_pMassOperator->reset(params->d_pMassOperatorParameters);
    }

  boost::shared_ptr<AMP::Operator::LinearOperator> pRhsOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearOperator>(d_pRhsOperator);
  AMP_INSIST(pRhsOperator.get()!=NULL, "ERROR: RhsOperator is not of type LinearOperator");
  
  boost::shared_ptr<AMP::LinearAlgebra::Matrix> pMatrix = pRhsOperator->getMatrix();

  if(d_bModifyRhsOperatorMatrix)
    {
      AMP_INSIST(pMatrix.get()!=NULL, "ERROR: NULL matrix pointer");
      setMatrix(pMatrix);
    }
  else
    {
      // if it's not okay to modify the rhs matrix then copy it over
      if(d_matrix.get()==NULL)
      {    
        AMP_INSIST(pMatrix.get()!=NULL, "ERROR: NULL matrix pointer");
        d_matrix = pMatrix->cloneMatrix();
        AMP_INSIST(d_matrix.get()!=NULL, "ERROR: NULL matrix pointer");
      }
      
      d_matrix->zero();
      d_matrix->makeConsistent();
      d_matrix->makeConsistent();

      d_matrix->axpy(1.0, *pMatrix);
    }

  if(!d_bAlgebraicComponent)
    {
      boost::shared_ptr<AMP::Operator::LinearOperator> pMassOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearOperator>(d_pMassOperator);
      AMP_INSIST(pMassOperator.get()!=NULL, "ERROR: MassOperator is not of type LinearOperator");
      
      boost::shared_ptr<AMP::LinearAlgebra::Matrix> pMassMatrix = pMassOperator->getMatrix();
      
      AMP_INSIST(pMassMatrix.get()!=NULL, "ERROR: NULL matrix pointer");
      // update the matrix to incorporate the contribution from the time derivative
      d_matrix->axpy(d_dScalingFactor, *pMassMatrix);
    }
        
  d_matrix->makeConsistent();
}

boost::shared_ptr<AMP::Operator::OperatorParameters>
LinearTimeOperator::getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& u)
{
  boost::shared_ptr<AMP::InputDatabase> timeOperator_db(new AMP::InputDatabase("LinearTimeOperatorDatabase"));
  timeOperator_db->putDouble("CurrentDt", d_dCurrentDt);
  timeOperator_db->putString("name", "LinearTimeOperator");
  timeOperator_db->putDouble("ScalingFactor", 1.0/d_dCurrentDt);

  boost::shared_ptr<AMP::TimeIntegrator::TimeOperatorParameters> timeOperatorParameters(new AMP::TimeIntegrator::TimeOperatorParameters(timeOperator_db));
  timeOperatorParameters->d_pRhsOperatorParameters = d_pRhsOperator->getJacobianParameters(u);
  timeOperatorParameters->d_pMassOperatorParameters = d_pMassOperator->getJacobianParameters(u);
  
  return timeOperatorParameters ;
}
  
}
}

