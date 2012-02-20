#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"
#include "materials/Material.h"

#include "ampmesh/MeshVariable.h"
#include "ampmesh/SiloIO.h"


#include "operators/MechanicsLinearFEOperator.h"
#include "operators/MechanicsNonlinearFEOperator.h"
#include "operators/DirichletVectorCorrection.h"
#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "../NonlinearKrylovAcceleratorParameters.h"
#include "../NonlinearKrylovAccelerator.h"

#include "../TrilinosMLSolver.h"



void myTest(AMP::UnitTest *ut, std::string exeName)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ( "brick" );

  AMP_INSIST(input_db->keyExists("NumberOfLoadingSteps"), "Key ''NumberOfLoadingSteps'' is missing!");
  int NumberOfLoadingSteps = input_db->getInteger("NumberOfLoadingSteps");

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  boost::shared_ptr<AMP::Database> nonlinOpDatabase = input_db->getDatabase("nonlinearMechanicsBVPOperator");
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinBvpOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter, 
          nonlinOpDatabase, elementPhysicsModel));
  (boost::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(nonlinBvpOperator->getVolumeOperator()))->init();

  boost::shared_ptr<AMP::Database> linOpDatabase = input_db->getDatabase("linearMechanicsBVPOperator");
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linBvpOperator =
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          linOpDatabase, elementPhysicsModel));

  AMP::LinearAlgebra::Variable::shared_ptr displacementVariable = boost::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
      nonlinBvpOperator->getVolumeOperator())->getInputVariable(AMP::Operator::Mechanics::DISPLACEMENT); 
  AMP::LinearAlgebra::Variable::shared_ptr residualVariable = nonlinBvpOperator->getOutputVariable();

  //For RHS (Point Forces)
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Database> load_db = input_db->getDatabase("Load_Boundary");
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletLoadVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          load_db, dummyModel));
  dirichletLoadVecOp->setVariable(residualVariable);

  //For Initial-Guess
  boost::shared_ptr<AMP::Database> disp_db = input_db->getDatabase("Displacement_Boundary");
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletDispInVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          disp_db, dummyModel));
  dirichletDispInVecOp->setVariable(displacementVariable);

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  AMP::LinearAlgebra::Vector::shared_ptr mechNlSolVec = meshAdapter->createVector( displacementVariable );
  AMP::LinearAlgebra::Vector::shared_ptr mechNlRhsVec = meshAdapter->createVector( residualVariable );
  AMP::LinearAlgebra::Vector::shared_ptr mechNlResVec = meshAdapter->createVector( residualVariable );
  AMP::LinearAlgebra::Vector::shared_ptr mechNlScaledRhsVec = meshAdapter->createVector( residualVariable );

#ifdef USE_SILO
  manager->registerVectorAsData ( mechNlSolVec , "Solution_Vector" );
  manager->registerVectorAsData ( mechNlResVec , "Residual_Vector" );
#endif

  //Initial guess for NL solver must satisfy the displacement boundary conditions
  mechNlSolVec->setToScalar(1.0);
  dirichletDispInVecOp->apply(nullVec, nullVec, mechNlSolVec, 1.0, 0.0);

  nonlinBvpOperator->apply(nullVec, mechNlSolVec, mechNlResVec, 1.0, 0.0);
  linBvpOperator->reset(nonlinBvpOperator->getJacobianParameters(mechNlSolVec));

  //Point forces
  mechNlRhsVec->setToScalar(0.0);
  dirichletLoadVecOp->apply(nullVec, nullVec, mechNlRhsVec, 1.0, 0.0);

  boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
  boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

  // ---- first initialize the preconditioner
  boost::shared_ptr<AMP::Database> pcSolver_db = linearSolver_db->getDatabase("Preconditioner"); 
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> pcSolverParams(new AMP::Solver::SolverStrategyParameters(pcSolver_db));
  pcSolverParams->d_pOperator = linBvpOperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(new AMP::Solver::TrilinosMLSolver(pcSolverParams));

  //HACK to prevent a double delete on Petsc Vec
  boost::shared_ptr<AMP::Solver::NonlinearKrylovAccelerator> nonlinearSolver;

  boost::shared_ptr<AMP::Solver::NonlinearKrylovAcceleratorParameters> nonlinearSolverParams(new
      AMP::Solver::NonlinearKrylovAcceleratorParameters(nonlinearSolver_db));

  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_pOperator = nonlinBvpOperator;
  nonlinearSolverParams->d_pPreconditioner = pcSolver;
  nonlinearSolverParams->d_pInitialGuess = mechNlSolVec;

  nonlinearSolver.reset(new AMP::Solver::NonlinearKrylovAccelerator(nonlinearSolverParams));

  nonlinearSolver->setZeroInitialGuess(false);

  for (int step=0;step<NumberOfLoadingSteps; step++)
  {
    AMP::pout << "########################################" << std::endl;
    AMP::pout << "The current loading step is " << (step+1) << std::endl;

    double scaleValue  = ((double)step+1.0)/NumberOfLoadingSteps;
    mechNlScaledRhsVec->scale(scaleValue, mechNlRhsVec);
    AMP::pout << "L2 Norm of RHS at loading step " << (step+1) << " is " << mechNlScaledRhsVec->L2Norm() << std::endl;

    nonlinBvpOperator->apply(mechNlScaledRhsVec, mechNlSolVec, mechNlResVec, 1.0, -1.0);
    double initialResidualNorm  = mechNlResVec->L2Norm();
    AMP::pout<<"Initial Residual Norm for loading step "<<(step+1)<<" is "<<initialResidualNorm<<std::endl;

    AMP::pout<<"Starting Nonlinear Solve..."<<std::endl;
    nonlinearSolver->solve(mechNlScaledRhsVec, mechNlSolVec);

    nonlinBvpOperator->apply(mechNlScaledRhsVec, mechNlSolVec, mechNlResVec, 1.0, -1.0);
    double finalResidualNorm  = mechNlResVec->L2Norm();
    AMP::pout<<"Final Residual Norm for loading step "<<(step+1)<<" is "<<finalResidualNorm<<std::endl;

    if( finalResidualNorm > (1.0e-10*initialResidualNorm) ) {
      ut.failure("Nonlinear solve for current loading step");
    } else {
      ut.passes("Nonlinear solve for current loading step");
    }

    boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
    boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperatorParameters> tmpParams(new
        AMP::Operator::MechanicsNonlinearFEOperatorParameters(tmp_db));
    (nonlinBvpOperator->getVolumeOperator())->reset(tmpParams);
    nonlinearSolver->setZeroInitialGuess(false);
  }

  double finalSolNorm = mechNlSolVec->L2Norm();
  AMP::pout<<"Final Solution Norm: "<<finalSolNorm<<std::endl;

#ifdef USE_SILO
  manager->writeFile<AMP::Mesh::SiloIO> ( "vonMisesExample" , 1 );
#endif

  ut.passes(exeName);

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    int node  = AMP::AMP_MPI::getRank();
    int nodes = AMP::AMP_MPI::getNodes();

    std::vector<std::string> exeNames;
    //  exeNames.push_back("testNonlinearKrylovAccelerator-LU-NonlinearMechanics-1-normal");
    exeNames.push_back("testNonlinearKrylovAccelerator-ML-NonlinearMechanics-1-normal");
    //  exeNames.push_back("testNonlinearKrylovAccelerator-LU-NonlinearMechanics-1-reduced");
    //  exeNames.push_back("testNonlinearKrylovAccelerator-ML-NonlinearMechanics-1-reduced");

    for(int i = 0; i < exeNames.size(); i++) {
        try {
            myTest(ut, exeNames[i]);
        } catch (std::exception &err) {
            AMP::pout << "ERROR: While testing "<<argv[0] 
                << err.what()
                << std::endl;
        } catch( ... ) {
            AMP::pout << "ERROR: While testing "<<argv[0] 
                << "An unknown exception was thrown."
                << std::endl;
        }
    } //end for i

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   


