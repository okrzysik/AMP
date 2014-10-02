#include "utils/InputDatabase.h"
#include "ColumnTimeOperator.h"
#include "TimeOperatorParameters.h"
#include "LinearTimeOperator.h"
#include "operators/ColumnOperatorParameters.h"

namespace AMP{
namespace TimeIntegrator{

ColumnTimeOperator::ColumnTimeOperator(boost::shared_ptr<AMP::Operator::OperatorParameters > in_params):ColumnOperator(in_params)
{
  
  boost::shared_ptr<TimeOperatorParameters> params = boost::dynamic_pointer_cast<TimeOperatorParameters>(in_params);
  boost::shared_ptr<AMP::Database> column_db = params->d_db;
  d_Mesh = params->d_Mesh;  
  d_pRhsOperator = boost::dynamic_pointer_cast<ColumnOperator>(params->d_pRhsOperator);
  AMP_INSIST(d_pRhsOperator.get()!=NULL, "Error: ColumnTimeOperator::ColumnTimeOperator() rhs operator must be a non-NULL column operator");

  d_pMassOperator = boost::dynamic_pointer_cast<ColumnOperator>(params->d_pMassOperator);
  AMP_INSIST(d_pRhsOperator.get()!=NULL, "Error: ColumnTimeOperator::ColumnTimeOperator() mass operator must be a non-NULL column operator");

  d_bCreateLinearTimeOperators = column_db->getBoolWithDefault("CreateLinearTimeOperators", true);

  // for now we only allow one algebraic component
   d_iAlgebraicComponent = column_db->getIntegerWithDefault("algebraicComponent", -1);

  if(d_bCreateLinearTimeOperators)
    {
      const int numberOfOperators = d_pRhsOperator->getNumberOfOperators();
      
      for( int i=0; i<numberOfOperators; i++)
    {
      boost::shared_ptr<AMP::InputDatabase> timeOperator_db(new AMP::InputDatabase("TimeOperatorDatabase"));
      // we assume for now that either all operators in the column operator are linear or all are nonlinear
      timeOperator_db->putDouble("CurrentDt", column_db->getDoubleWithDefault("CurrentDt", 1.0e-08));
      timeOperator_db->putDouble("CurrentTime", column_db->getDoubleWithDefault("CurrentTime", 0.0));
      timeOperator_db->putString("name", "TimeOperator");
      timeOperator_db->putBool("bLinearMassOperator", column_db->getBoolWithDefault("bLinearMassOperator", true));
      timeOperator_db->putBool("bLinearRhsOperator", column_db->getBoolWithDefault("bLinearRhsOperator", false));
      timeOperator_db->putDouble("ScalingFactor", column_db->getDoubleWithDefault("ScalingFactor", 1.0e6));
      boost::shared_ptr<AMP::TimeIntegrator::TimeOperatorParameters> timeOperatorParameters(new AMP::TimeIntegrator::TimeOperatorParameters(timeOperator_db));
      timeOperatorParameters->d_pRhsOperator = d_pRhsOperator->getOperator(i);

      // if there are algebraic components set the mass operator to NULL
      if( i != d_iAlgebraicComponent)
        {
          timeOperatorParameters->d_pMassOperator = d_pMassOperator->getOperator(i);
        }
      else
        {
          timeOperator_db->putBool("bAlgebraicComponent", true);
        }

      boost::shared_ptr< AMP::TimeIntegrator::LinearTimeOperator > op(new AMP::TimeIntegrator::LinearTimeOperator(timeOperatorParameters));

      d_Operators.push_back( op );

    }
    }
  else
    {
        AMP::pout << "Error: ColumnTimeOperator::ColumnTimeOperator() currently implemented for column rhs and mass operators where all components are LinearOperators" << std::endl;
    }

  reset(in_params);

  d_pSourceTerm   = params->d_pSourceTerm;
  d_pPreviousTimeSolution = params->d_pPreviousTimeSolution;

}
  
ColumnTimeOperator::~ColumnTimeOperator()
{
}

void
ColumnTimeOperator::reset(const boost::shared_ptr<AMP::Operator::OperatorParameters>& in_params)
{
  boost::shared_ptr<TimeOperatorParameters> params = boost::dynamic_pointer_cast<TimeOperatorParameters>(in_params);
  boost::shared_ptr<AMP::Database> column_db = params->d_db;
  boost::shared_ptr<AMP::Operator::ColumnOperatorParameters> pRhsParameters =  boost::dynamic_pointer_cast<AMP::Operator::ColumnOperatorParameters>(params->d_pRhsOperatorParameters);
  boost::shared_ptr<AMP::Operator::ColumnOperatorParameters> pMassParameters =  boost::dynamic_pointer_cast<AMP::Operator::ColumnOperatorParameters>(params->d_pMassOperatorParameters);
  
  AMP_INSIST(params.get()!=NULL, "Error: NULL TimeOperatorParameters object");
  
  getFromInput(params->d_db);

  const int numberOfOperators = d_pRhsOperator->getNumberOfOperators();
  
  for( int i=0; i<numberOfOperators; i++)
    {
      boost::shared_ptr<AMP::InputDatabase> timeOperator_db(new AMP::InputDatabase("TimeOperatorDatabase"));
      // we assume for now that either all operators in the column operator are linear or all are nonlinear
      timeOperator_db->putDouble("CurrentDt", column_db->getDoubleWithDefault("CurrentDt", 1.0e-08));
      timeOperator_db->putDouble("CurrentTime", column_db->getDoubleWithDefault("CurrentTime", 0.0));
      timeOperator_db->putString("name", "TimeOperator");
      timeOperator_db->putBool("bLinearMassOperator", column_db->getBoolWithDefault("bLinearMassOperator", true));
      timeOperator_db->putBool("bLinearRhsOperator", column_db->getBoolWithDefault("bLinearRhsOperator", false));
      timeOperator_db->putDouble("ScalingFactor", column_db->getDoubleWithDefault("ScalingFactor", 1.0e6));
      boost::shared_ptr<AMP::TimeIntegrator::TimeOperatorParameters> timeOperatorParameters(new AMP::TimeIntegrator::TimeOperatorParameters(timeOperator_db));
      if(pRhsParameters.get()!=NULL)
    {
      timeOperatorParameters->d_pRhsOperatorParameters = (pRhsParameters->d_OperatorParameters)[i];
    }

      // if there are algebraic components set the mass operator to NULL
      if( i != d_iAlgebraicComponent )
    {
      if(pMassParameters.get()!=NULL)
        {
          timeOperatorParameters->d_pMassOperatorParameters = (pMassParameters->d_OperatorParameters)[i];
        }
    }
      else
    {
      timeOperator_db->putBool("bAlgebraicComponent", true);
    }
      
      d_Operators[i]->reset(timeOperatorParameters);
    }
}
  
void
ColumnTimeOperator::getFromInput(const boost::shared_ptr<AMP::Database> &db)
{
  AMP_INSIST(db->keyExists("CurrentDt"), "key CurrentDt missing in input");

  d_dCurrentDt = db->getDouble("CurrentDt");
  
}
  
void
ColumnTimeOperator::apply(AMP::LinearAlgebra::Vector::const_shared_ptr /* f */, AMP::LinearAlgebra::Vector::const_shared_ptr /* u */, 
              AMP::LinearAlgebra::Vector::shared_ptr /* r */, const double /* a */, const double /* b */)
{
  abort();
}

void
ColumnTimeOperator::append(boost::shared_ptr< Operator > /* op */)
{
  AMP::pout << "Error: ColumnTimeOperator::append(): this routine is disabled for ColumnTimeOperators" << std::endl;
}
  
}
}

