#include "solvers/trilinos/nox/TrilinosNOXSolver.h"
#include "utils/ProfilerApp.h"

#include "vectors/trilinos/thyra/ThyraVector.h"
#include "vectors/trilinos/thyra/ThyraVectorWrapper.h"
#include "solvers/trilinos/thyra/TrilinosThyraModelEvaluator.h"


// Trilinos includes
#include "NOX_Thyra.H"
#include "NOX_Thyra_Group.H"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "NOX_StatusTest_Combo.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_StatusTest_FiniteValue.H"
#include "NOX_Solver_Factory.H"
#include "BelosTypes.hpp"
#include "NOX_Thyra_MatrixFreeJacobianOperator.hpp"
#include "NOX_MatrixFree_ModelEvaluatorDecorator.hpp"


namespace AMP {
namespace Solver {


/****************************************************************
*  Constructors                                                 *
****************************************************************/
TrilinosNOXSolver::TrilinosNOXSolver():
    SolverStrategy()
{
    
}
TrilinosNOXSolver::TrilinosNOXSolver(boost::shared_ptr<TrilinosNOXSolverParameters> parameters):
    SolverStrategy(parameters)
{
    TrilinosNOXSolver();
    initialize(parameters);
}
void TrilinosNOXSolver::reset(boost::shared_ptr<SolverStrategyParameters> parameters) 
{
    initialize(parameters);
}
TrilinosNOXSolver::~TrilinosNOXSolver()
{
}


/****************************************************************
*  Initialize                                                   *
****************************************************************/
void TrilinosNOXSolver::initialize( boost::shared_ptr<SolverStrategyParameters> parameters )
{
    // Copy the parameters
    boost::shared_ptr<TrilinosNOXSolverParameters> params = 
        boost::dynamic_pointer_cast<TrilinosNOXSolverParameters>( parameters );
    AMP_ASSERT(params.get()!=NULL);
    if ( params->d_pInitialGuess.get()!=NULL )
        d_initialGuess = params->d_pInitialGuess;
    AMP_ASSERT(d_initialGuess!=NULL);
    boost::shared_ptr<AMP::Database> nonlinear_db = parameters->d_db;
    boost::shared_ptr<AMP::Database> linear_db = nonlinear_db->getDatabase("LinearSolver");
    AMP_ASSERT(linear_db!=NULL);
    // Create the default OStream
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    // Create a model evaluator
    boost::shared_ptr<TrilinosThyraModelEvaluatorParameters> modelParams( new TrilinosThyraModelEvaluatorParameters );
    modelParams->d_nonlinearOp = d_pOperator;
    modelParams->d_linearOp = params->d_pLinearOperator;
    modelParams->d_icVec = d_initialGuess;
    //modelParams->d_preconditioner.reset();
    modelParams->d_preconditioner = params->d_preconditioner;
    d_thyraModel = Teuchos::RCP<TrilinosThyraModelEvaluator>( new TrilinosThyraModelEvaluator(modelParams) );
    // Create the Preconditioner operator
    d_precOp = d_thyraModel->create_W_prec();
    // Create the linear solver factory
    ::Stratimikos::DefaultLinearSolverBuilder builder;
    Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
    std::string linearSolverType = linear_db->getString("linearSolverType");
    std::string linearSolver     = linear_db->getString("linearSolver");
    p->set("Linear Solver Type",linearSolverType);
    p->set("Preconditioner Type","None");
    p->sublist("Linear Solver Types").sublist(linearSolverType).set("Solver Type",linearSolver);
    Teuchos::ParameterList& linearSolverParams = p->sublist("Linear Solver Types").sublist(linearSolverType);
linearSolverParams.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Maximum Iterations",100);
    if ( linear_db->getIntegerWithDefault("print_info_level",0) >= 2 ) {
        linearSolverParams.sublist("Solver Types").sublist(linearSolver).set("Output Frequency",1);
        linearSolverParams.sublist("Solver Types").sublist(linearSolver).set("Verbosity",10);
        linearSolverParams.sublist("VerboseObject").set("Verbosity Level","extreme");
        if ( linearSolverType == "Belos" ) {
            linearSolverParams.sublist("Solver Types").sublist(linearSolver).set("Verbosity",
                Belos::Warnings+Belos::IterationDetails+Belos::OrthoDetails+
                Belos::FinalSummary+Belos::Debug+Belos::StatusTestDetails);
        }
    }
    //Teuchos::ParameterList::writeParameterListToXmlFile(*p,xmlFileName);
    builder.setParameterList(p);
    d_lowsFactory = builder.createLinearSolveStrategy("");
    d_lowsFactory->initializeVerboseObjectBase();
    d_thyraModel->set_W_factory(d_lowsFactory);
    // Create the convergence tests (these will need to be on the input database)
    Teuchos::RCP<NOX::StatusTest::NormF> absresid =
        Teuchos::rcp(new NOX::StatusTest::NormF(d_dMaxError));
    Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
        Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2,d_dMaxError));
    Teuchos::RCP<NOX::StatusTest::Combo> converged =
        Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
    converged->addStatusTest(absresid);
    converged->addStatusTest(wrms);
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
        Teuchos::rcp(new NOX::StatusTest::MaxIters(d_iMaxIterations));
    Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
        Teuchos::rcp(new NOX::StatusTest::FiniteValue);
    d_status = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    d_status->addStatusTest(fv);
    d_status->addStatusTest(converged);
    d_status->addStatusTest(maxiters);
    // Create nox parameter list
    d_nlParams = Teuchos::rcp(new Teuchos::ParameterList);
    d_nlParams->set("Nonlinear Solver", "Line Search Based");
    // Set the printing parameters in the "Printing" sublist
    Teuchos::ParameterList& printParams = d_nlParams->sublist("Printing");
    printParams.set("Output Precision", 3);
    printParams.set("Output Processor", 0);
    if ( d_iDebugPrintInfoLevel >= 2 ) {
        printParams.set("Output Information", 
			     NOX::Utils::OuterIteration + 
			     NOX::Utils::OuterIterationStatusTest + 
			     NOX::Utils::InnerIteration +
			     NOX::Utils::LinearSolverDetails +
			     NOX::Utils::Parameters + 
			     NOX::Utils::Details + 
			     NOX::Utils::Warning +
                             NOX::Utils::Debug +
			     NOX::Utils::TestDetails +
			     NOX::Utils::Error);
    }
}


/****************************************************************
*  Solve                                                        *
****************************************************************/
void TrilinosNOXSolver::solve( boost::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                  boost::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE_START("solve");
    // Get thyra vectors
    boost::shared_ptr<AMP::LinearAlgebra::ThyraVector> initial = 
        boost::dynamic_pointer_cast<AMP::LinearAlgebra::ThyraVector>(
        AMP::LinearAlgebra::ThyraVector::view( d_initialGuess ) );
    boost::shared_ptr<AMP::LinearAlgebra::ThyraVector> U = 
        boost::dynamic_pointer_cast<AMP::LinearAlgebra::ThyraVector>(
        AMP::LinearAlgebra::ThyraVector::view( u ) );
    boost::shared_ptr<const AMP::LinearAlgebra::ThyraVector> F = 
        boost::dynamic_pointer_cast<const AMP::LinearAlgebra::ThyraVector>(
        AMP::LinearAlgebra::ThyraVector::constView( f ) );
    // Set the rhs for the thyra model
    d_thyraModel->setRhs( f );
    // Create the JFNK operator
    Teuchos::ParameterList printParams;
    Teuchos::RCP<Teuchos::ParameterList> jfnkParams = Teuchos::parameterList();
    jfnkParams->set("Difference Type","Forward");
    jfnkParams->set("Perturbation Algorithm","KSP NOX 2001");
    jfnkParams->set("lambda",1.0e-4);
    Teuchos::RCP<NOX::Thyra::MatrixFreeJacobianOperator<double> > jfnkOp(
        new NOX::Thyra::MatrixFreeJacobianOperator<double>(printParams) );
    jfnkOp->setParameterList(jfnkParams);
    jfnkParams->print(std::cout);
    // Create the NOX::Thyra::Group
    //Teuchos::RCP<NOX::Thyra::Group> nox_group( new NOX::Thyra::Group( initial->getVec(), d_thyraModel ) );
    Teuchos::RCP< ::Thyra::ModelEvaluator<double> > thyraModel = 
        Teuchos::rcp(new NOX::MatrixFreeModelEvaluatorDecorator<double>(d_thyraModel));
    Teuchos::RCP<NOX::Thyra::Group> nox_group( new NOX::Thyra::Group( initial->getVec(), thyraModel, jfnkOp, d_lowsFactory, d_precOp, Teuchos::null));
    nox_group->setX(U->getVec());
    nox_group->computeF();
    // VERY IMPORTANT!!!  jfnk object needs base evaluation objects.
    // This creates a circular dependency, so use a weak pointer.
    jfnkOp->setBaseEvaluationToNOXGroup(nox_group.create_weak());
    // Create the solver
    d_solver = NOX::Solver::buildSolver(nox_group, d_status, d_nlParams);
    // Solve
    NOX::StatusTest::StatusType solvStatus = d_solver->solve();
    if ( solvStatus != NOX::StatusTest::Converged )
        AMP_ERROR("Failed to solve");
    // Copy the solution back to u
    const NOX::Thyra::Vector* tmp = dynamic_cast<const NOX::Thyra::Vector*>(&(nox_group->getX()));
    const AMP::LinearAlgebra::ThyraVectorWrapper* thyraVec = 
        dynamic_cast<const AMP::LinearAlgebra::ThyraVectorWrapper*>(&(tmp->getThyraVector()));
    AMP_ASSERT(thyraVec!=NULL);
    AMP_ASSERT(thyraVec->numVecs()==1);
    u->copyVector(thyraVec->getVec(0));
    PROFILE_STOP("solve");
}




}
}

