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


namespace AMP {
namespace Solver {


/****************************************************************
*  Constructors                                                 *
****************************************************************/
TrilinosNOXSolver::TrilinosNOXSolver()
{
    
}
TrilinosNOXSolver::TrilinosNOXSolver(boost::shared_ptr<TrilinosNOXSolverParameters> parameters):SolverStrategy(parameters)
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
    // Create the default OStream
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    //Teuchos::RCP<Teuchos::FancyOStream> defaultOStream( new Teuchos::FancyOStream(Teuchos::rcpFromRef(std::cout)) );
    //Teuchos::VerboseObjectBase::setDefaultOStream(defaultOStream);
    // Create a model evaluator
    boost::shared_ptr<TrilinosThyraModelEvaluatorParameters> modelParams( new TrilinosThyraModelEvaluatorParameters );
    modelParams->d_nonlinearOp = d_pOperator;
    modelParams->d_linearOp = params->d_pLinearOperator;
    modelParams->d_icVec = d_initialGuess;
    d_thyraModel = Teuchos::RCP<TrilinosThyraModelEvaluator>( new TrilinosThyraModelEvaluator(modelParams) );
    // Create the linear solver factory
    ::Stratimikos::DefaultLinearSolverBuilder builder;
    Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
    p->set("Linear Solver Type", "Belos");
    p->set("Preconditioner Type", "None");
    //p->set("Linear Solver Type", "AztecOO");
    //p->set("Preconditioner Type", "Ifpack");
    //p->set("Enable Delayed Solver Construction", true);
    Teuchos::ParameterList& linearSolverParams = p->sublist("Linear Solver Types");
    Teuchos::ParameterList& belosParams = linearSolverParams.sublist("Belos");
    std::string belosSolverType = "Block CG";
    belosParams.set("Solver Type",belosSolverType);
    Teuchos::ParameterList& belosSolverTypeParams = belosParams.sublist("Solver Types");
    Teuchos::ParameterList& belosSolverParams = belosSolverTypeParams.sublist(belosSolverType);
    if ( d_iDebugPrintInfoLevel >= 2 ) {
        belosSolverParams.set("Verbosity",Belos::Warnings+Belos::IterationDetails+
            Belos::OrthoDetails+Belos::FinalSummary+Belos::Debug+Belos::StatusTestDetails);
    }
    builder.setParameterList(p);
    Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> > 
        lowsFactory = builder.createLinearSolveStrategy("");
    lowsFactory->initializeVerboseObjectBase();
    d_thyraModel->set_W_factory(lowsFactory);
    // Create the convergence tests (these will need to be on the input database)
    Teuchos::RCP<NOX::StatusTest::NormF> absresid =
        Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
    Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
        Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
    Teuchos::RCP<NOX::StatusTest::Combo> converged =
        Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
    converged->addStatusTest(absresid);
    converged->addStatusTest(wrms);
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
        Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
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
    //printParams.set("MyPID",0); 
    printParams.set("Output Precision", 3);
    printParams.set("Output Processor", 0);
    //printParams.set("Output Information",NOX::Utils::Error+NOX::Utils::TestDetails);
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
void TrilinosNOXSolver::solve( boost::shared_ptr<AMP::LinearAlgebra::Vector> f,
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
    boost::shared_ptr<AMP::LinearAlgebra::ThyraVector> F = 
        boost::dynamic_pointer_cast<AMP::LinearAlgebra::ThyraVector>(
        AMP::LinearAlgebra::ThyraVector::view( f ) );
    // Set the rhs for the thyra model
    d_thyraModel->setRhs( f );
    // Create the NOX::Thyra::Group
    Teuchos::RCP<NOX::Thyra::Group> nox_group( new NOX::Thyra::Group( initial->getVec(), d_thyraModel ) );
    nox_group->setX(U->getVec());
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

