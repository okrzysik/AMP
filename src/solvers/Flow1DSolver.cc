#include "utils/Utilities.h"
#include "FlowFrapconJacobian.h"
#include "Flow1DSolver.h"
#include "MultiVector.h"

namespace AMP {
namespace Solver {

Flow1DSolver::Flow1DSolver(boost::shared_ptr<SolverStrategyParameters> parameters):SolverStrategy(parameters)
{
  
  assert(parameters.get()!=NULL);
  
  initialize(parameters);

}

Flow1DSolver::~Flow1DSolver()
{

}

void
Flow1DSolver::setInitialGuess( boost::shared_ptr<AMP::LinearAlgebra::Vector>   )
{
  
}

void
Flow1DSolver::initialize(boost::shared_ptr<SolverStrategyParameters> const parameters)
{
  getFromInput(parameters->d_db);

  if(d_pOperator.get()!=NULL)
    {
      boost::shared_ptr<AMP::Operator::FlowFrapconJacobian> Operator = boost::dynamic_pointer_cast<AMP::Operator::FlowFrapconJacobian>(d_pOperator);

      assert(Operator.get() != NULL);

      AMP_ERROR("Flow1DSolver is not converted yet");
      //d_inpVariable = Operator->d_inpVariable;
      //d_outVariable = Operator->d_outVariable;
      d_K =  Operator->d_K;
      d_De = Operator->d_De;
      d_G = Operator->d_G;
      d_Re = Operator->d_Re;
      d_Pr = Operator->d_Pr;
      d_Cp = Operator->getCp();
      d_dCp = Operator->getdCp();
      zPoints = Operator->getZLocations();
      d_cladVec = Operator->d_cladVec;

    }
}

void
Flow1DSolver::reset(boost::shared_ptr<SolverStrategyParameters> )
{

  if(d_pOperator.get()!=NULL)
    {
    }
}

void
Flow1DSolver::resetOperator(const boost::shared_ptr<AMP::Operator::OperatorParameters> params)
{
  if(d_pOperator.get()!=NULL)
    {
      d_pOperator->reset(params);
    }
}

  void
Flow1DSolver::solve(boost::shared_ptr<AMP::LinearAlgebra::Vector>  f,
		          boost::shared_ptr<AMP::LinearAlgebra::Vector>  u)
{
  
  AMP::LinearAlgebra::Vector::shared_ptr flowInputVec = u->subsetVectorForVariable(d_inpVariable);
  AMP::LinearAlgebra::Vector::shared_ptr flowRhsVec   = f->subsetVectorForVariable(d_inpVariable);

  d_numpoints   = zPoints.size();

  for( int i=1; i<d_numpoints; i++) {

    double cur_node, next_node;

    cur_node  = zPoints[i-1];
    next_node = zPoints[i];

    // double Heff, he_z, T_b_i, T_b_im1, T_c_i;
    // double T_b =0, Tin, frhs;
    double Heff, he_z, T_b_i, T_b_im1;
    double frhs;

    Heff = (0.023*d_K/d_De)*pow(d_Re,0.8)*pow(d_Pr,0.4);
    he_z = next_node - cur_node;

    // T_c_i   = d_cladVec->getValueByLocalID(i);
    T_b_im1 = flowInputVec->getValueByLocalID(i-1);
    T_b_i   = flowInputVec->getValueByLocalID(i);
    frhs    = flowRhsVec->getValueByLocalID(i);

    T_b_i = ( T_b_im1 + frhs ) * (1./(1.0 + ((4*Heff*he_z)/(d_Cp*d_G*d_De))));

    flowInputVec->setValueByLocalID(i, T_b_i);

  }//end for i
}

}
}

