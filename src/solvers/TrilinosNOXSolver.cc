#include "solvers/TrilinosNOXSolver.h"
#include "utils/ProfilerApp.h"

#include "vectors/trilinos/ThyraVector.h"
#include "solvers/TrilinosThyraModelEvaluator.h"


// Trilinos includes
#include "NOX_Thyra.H"
#include "NOX_Thyra_Group.H"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"


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
    Teuchos::RCP<TrilinosThyraModelEvaluator> thyraModel( new TrilinosThyraModelEvaluator() );
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
}


/****************************************************************
*  Solve                                                        *
****************************************************************/
void TrilinosNOXSolver::solve( boost::shared_ptr<AMP::LinearAlgebra::Vector> f,
                  boost::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE_START("solve");
    
    
    
    PROFILE_STOP("solve");
}




}
}

