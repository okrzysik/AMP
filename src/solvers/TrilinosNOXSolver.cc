#include "solvers/TrilinosNOXSolver.h"
#include "utils/ProfilerApp.h"

#include "vectors/trilinos/ThyraVector.h"
#include "solvers/TrilinosThyraModelEvaluator.h"


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


namespace AMP {
namespace Solver {


/****************************************************************
*  Constructors                                                 *
****************************************************************/
TrilinosNOXSolver::TrilinosNOXSolver()
{
    
}
TrilinosNOXSolver::TrilinosNOXSolver(boost::shared_ptr<TrilinosNOXSolverParameters> parameters)
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
    // Create a model evaluator
    boost::shared_ptr<TrilinosThyraModelEvaluatorParameters> modelParams( new TrilinosThyraModelEvaluatorParameters );
    modelParams->d_dofs = d_initialGuess->getDOFManager();
    Teuchos::RCP<TrilinosThyraModelEvaluator> thyraModel( new TrilinosThyraModelEvaluator(modelParams) );
    // Create the linear solver factory
    ::Stratimikos::DefaultLinearSolverBuilder builder;
    Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
    p->set("Linear Solver Type", "AztecOO");
    p->set("Preconditioner Type", "Ifpack");
    //p->set("Enable Delayed Solver Construction", true);
    builder.setParameterList(p);
    Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> > 
        lowsFactory = builder.createLinearSolveStrategy("");
    thyraModel->set_W_factory(lowsFactory);
    // Get the initial guess
    boost::shared_ptr<AMP::LinearAlgebra::ThyraVector> thyraVec = 
        boost::dynamic_pointer_cast<AMP::LinearAlgebra::ThyraVector>(
        AMP::LinearAlgebra::ThyraVector::view( d_initialGuess ) );
    Teuchos::RCP< ::Thyra::VectorBase<double> > initial_guess = thyraVec->getVec();
    // Create the NOX::Thyra::Group
    Teuchos::RCP<NOX::Thyra::Group> nox_group( new NOX::Thyra::Group( *initial_guess, thyraModel ) );
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
    Teuchos::RCP<NOX::StatusTest::Combo> combo =
        Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(fv);
    combo->addStatusTest(converged);
    combo->addStatusTest(maxiters);
    // Create nox parameter list
    Teuchos::RCP<Teuchos::ParameterList> nlParams =
        Teuchos::rcp(new Teuchos::ParameterList);
    nlParams->set("Nonlinear Solver", "Line Search Based");
    // Set the printing parameters in the "Printing" sublist
    Teuchos::ParameterList& printParams = nlParams->sublist("Printing");
    //printParams.set("MyPID",0); 
    printParams.set("Output Precision", 3);
    printParams.set("Output Processor", 0);
    printParams.set("Output Information",NOX::Utils::Error+NOX::Utils::TestDetails);
    // Create the solver
    d_solver = NOX::Solver::buildSolver(nox_group, combo, nlParams);

}


/****************************************************************
*  Solve                                                        *
****************************************************************/
void TrilinosNOXSolver::solve( boost::shared_ptr<AMP::LinearAlgebra::Vector> f,
                  boost::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE_START("solve");
    // Create the NOX::Solver
    NOX::StatusTest::StatusType solvStatus = d_solver->solve();
    if (solvStatus != NOX::StatusTest::Converged)
        AMP_ERROR("Failed to solve");
    PROFILE_STOP("solve");
}




}
}

