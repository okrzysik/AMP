#include "utils/Utilities.h"
#include "FlowFrapconJacobian.h"
#include "CoupledFlow1DSolver.h"
#include "Flow1DSolver.h"
#include "MultiVector.h"
#include "utils/InputDatabase.h"


namespace AMP {
namespace Solver {

CoupledFlow1DSolver::CoupledFlow1DSolver(boost::shared_ptr<SolverStrategyParameters> parameters):SolverStrategy(parameters)
{
AMP_ERROR("CoupledFlow1DSolver is not converted yet");
/*
  assert(parameters.get()!=NULL);

    boost::shared_ptr<CoupledFlow1DSolverParameters> params =
        boost::dynamic_pointer_cast<CoupledFlow1DSolverParameters>(parameters);

        std::string flowOutVar =( (params->d_pOperator)->getOutputVariable())->getName() ;

        d_pOperator = boost::dynamic_pointer_cast<AMP::Operator::CoupledFlowFrapconOperator>(params->d_pOperator);

        d_flow1DSolver = boost::dynamic_pointer_cast<Flow1DSolver>(params->d_flow1DSolver);
        std::string flowInpVar = (d_flow1DSolver->getInputVariable())->getName(); 

 	boost::shared_ptr<AMP::InputDatabase> tmp_db1 (new AMP::InputDatabase("Dummy"));
        tmp_db1->putInteger("BoundaryId",4);
        tmp_db1->putString("InputVariable", flowOutVar );
        tmp_db1->putString("OutputVariable", flowInpVar );
	boost::shared_ptr<AMP::Operator::MapOperatorParameters> mapflowInternal3to1Params (new AMP::Operator::MapOperatorParameters( tmp_db1 ));
	mapflowInternal3to1Params->d_MeshAdapter = (d_flow1DSolver->getOperator())->getMeshAdapter();
	d_flowInternal3to1.reset(new AMP::Operator::Map3Dto1D(mapflowInternal3to1Params) );

	boost::shared_ptr<AMP::InputDatabase> tmp_db2 (new AMP::InputDatabase("Dummy"));
        tmp_db2->putInteger("BoundaryId",4);
        tmp_db2->putString("InputVariable", flowInpVar );
        tmp_db2->putString("OutputVariable", flowOutVar );
	boost::shared_ptr<AMP::Operator::MapOperatorParameters> mapflowInternal1to3Params (new AMP::Operator::MapOperatorParameters( tmp_db2 ));
	mapflowInternal1to3Params->d_MapAdapter = (d_flow1DSolver->getOperator())->getMeshAdapter();
	d_flowInternal1to3.reset(new AMP::Operator::Map1Dto3D(mapflowInternal1to3Params) );

	(boost::dynamic_pointer_cast<AMP::Operator::Map3Dto1D> (d_flowInternal3to1) )->setZLocations( (boost::dynamic_pointer_cast<AMP::Operator::Map1Dto3D> (d_flowInternal1to3) )->getZLocations());

	int d_numpoints = (boost::dynamic_pointer_cast<AMP::Operator::Map1Dto3D> (d_flowInternal1to3) )->getNumZlocations();
        d_SimpleVariable.reset(new AMP::LinearAlgebra::Variable(flowInpVar));

	d_flowInput = AMP::LinearAlgebra::SimpleVector::create( d_numpoints, d_SimpleVariable ); 
	d_flowOutput = AMP::LinearAlgebra::SimpleVector::create( d_numpoints, d_SimpleVariable ); 
*/
}

CoupledFlow1DSolver::~CoupledFlow1DSolver()
{

}

void
CoupledFlow1DSolver::setInitialGuess( boost::shared_ptr<AMP::LinearAlgebra::Vector>   )
{
  
}

void
CoupledFlow1DSolver::reset(boost::shared_ptr<SolverStrategyParameters> )
{

  if(d_pOperator.get()!=NULL)
    {
    }
}

void
CoupledFlow1DSolver::resetOperator(const boost::shared_ptr<AMP::Operator::OperatorParameters> params)
{
  if(d_pOperator.get()!=NULL)
    {
      d_pOperator->reset(params);
    }
}

  void
CoupledFlow1DSolver::solve(boost::shared_ptr<AMP::LinearAlgebra::Vector>  f,
		          boost::shared_ptr<AMP::LinearAlgebra::Vector>  u)
{
AMP_ERROR("CoupledFlow1DSolver is not converted yet");
/*
	AMP::LinearAlgebra::Vector::shared_ptr   nullVec;
  
        d_inpVariable = d_flow1DSolver->getOperator()->getInputVariable();
        d_outVariable = d_flowInternal1to3->getOutputVariable();

        d_Sol = u->subsetVectorForVariable(d_outVariable);
        d_Rhs = f->subsetVectorForVariable(d_outVariable);

        (boost::dynamic_pointer_cast<AMP::Operator::Map3Dto1D> (d_flowInternal3to1))->setVector(d_flowInput);
        (boost::dynamic_pointer_cast<AMP::Operator::Map1Dto3D> (d_flowInternal1to3))->setVector(d_Sol );

        d_flowInternal3to1->apply(nullVec, d_Rhs, nullVec, -1, 1);
        d_flow1DSolver->solve(d_flowInput, d_flowOutput);
        d_flowInternal1to3->apply(nullVec, d_flowOutput , nullVec, -1, 1);
*/
}

}
}

