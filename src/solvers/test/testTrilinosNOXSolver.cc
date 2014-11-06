// This tests checks the creation of a TrilinosNOXSolver
#include <iostream>
#include <string>
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"
#include "vectors/NullVector.h"
#include "vectors/SimpleVector.h"
#include "vectors/MultiVector.h"
#include "operators/NullOperator.h"
#include "operators/IdentityOperator.h"
#include "solvers/trilinos/nox/TrilinosNOXSolver.h"



void myTest(AMP::UnitTest *ut, std::string exeName)
{
    std::string input_file = "input_" + exeName;
    std::string log_file = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero(log_file);
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    AMP::AMP_MPI solverComm = globalComm.dup();     // Create a unique solver comm to test proper cleanup

    AMP::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
    AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
    input_db->printClassData(AMP::plog);

    // Create the solution and function variables
    AMP::LinearAlgebra::Variable::shared_ptr var(new AMP::LinearAlgebra::Variable("x"));
    AMP::LinearAlgebra::Vector::shared_ptr u = AMP::LinearAlgebra::SimpleVector::create(25,var,solverComm);
    AMP::LinearAlgebra::Vector::shared_ptr f = u->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr icVec = u->cloneVector();

    // Create the operator
    AMP::shared_ptr<AMP::Operator::IdentityOperator> op(new AMP::Operator::IdentityOperator());
    op->setInputVariable(var);
    op->setOutputVariable(var);

    // Get the databases for the nonlinear and linear solvers
    AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
    //AMP::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

    // initialize the nonlinear solver parameters
    AMP::shared_ptr<AMP::Solver::TrilinosNOXSolverParameters> nonlinearSolverParams(new
       AMP::Solver::TrilinosNOXSolverParameters(nonlinearSolver_db));
    nonlinearSolverParams->d_comm = solverComm;
    nonlinearSolverParams->d_pInitialGuess = icVec;
    nonlinearSolverParams->d_pOperator = op;
    nonlinearSolverParams->d_pLinearOperator = op;

    // Create the nonlinear solver
    AMP::shared_ptr<AMP::Solver::TrilinosNOXSolver> nonlinearSolver(new AMP::Solver::TrilinosNOXSolver(nonlinearSolverParams));
    ut->passes("TrilinosNOXSolver created");

    // Call solve with a simple vector
    u->setRandomValues();
    f->setRandomValues();
    nonlinearSolver->solve(f,u);
    ut->passes("TrilinosNOXSolver solve called with simple vector");
    AMP::LinearAlgebra::Vector::shared_ptr x = u->cloneVector();
    x->subtract(u,f);
    double error  = x->L2Norm()/std::max(f->L2Norm(),1.0);
    if ( fabs(error)<1e-8 )
        ut->passes("Solve with simple vector passed");
    else
        ut->failure("Solve with simple vector failed");
    
    
    // Call solve with a multivector (there can be bugs when solve is called with a single vector and then a multivector)
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> mu = AMP::LinearAlgebra::MultiVector::create("multivector",solverComm);
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> mf = AMP::LinearAlgebra::MultiVector::create("multivector",solverComm);
    mu->addVector(u);
    mf->addVector(f);
    mu->setRandomValues();
    mf->zero();
    nonlinearSolver->solve(mf,mu);
    ut->passes("TrilinosNOXSolver solve called with multivector");

}



int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    myTest( &ut, "testTrilinosNOXSolver" );
   
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   

