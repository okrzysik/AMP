
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/PIO.h"

#include <iostream>
#include <string>

#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "utils/Writer.h"

#include "operators/mechanics/ThermalVonMisesMatModel.h"
#include "operators/mechanics/MechanicsNonlinearElement.h"
#include "operators/mechanics/MechanicsLinearElement.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"

#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscSNESSolverParameters.h"
#include "solvers/petsc/PetscSNESSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"

void myTest(AMP::UnitTest *ut)
{
  std::string exeName("testPetscSNESSolver-NonlinearMechanics-2");
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

#ifdef USE_EXT_SILO
  // Create the silo writer and register the data
  AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
#endif

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);

  AMP_INSIST(input_db->keyExists("NumberOfLoadingSteps"), "Key ''NumberOfLoadingSteps'' is missing!");
  int NumberOfLoadingSteps = input_db->getInteger("NumberOfLoadingSteps");

  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinBvpOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,"nonlinearMechanicsBVPOperator",input_db));
  boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> nonlinearMechanicsVolumeOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(nonlinBvpOperator->getVolumeOperator());
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel = nonlinearMechanicsVolumeOperator->getMaterialModel();

  AMP::LinearAlgebra::Variable::shared_ptr temperatureVariable(new AMP::LinearAlgebra::Variable("temp")); 
  AMP::LinearAlgebra::Variable::shared_ptr displacementVariable = nonlinBvpOperator->getOutputVariable();

  AMP::Discretization::DOFManager::shared_ptr tempDofMap = AMP::Discretization::simpleDOFManager::create(
      meshAdapter, AMP::Mesh::Vertex, 1, 1, true); 

  AMP::LinearAlgebra::Vector::shared_ptr initTempVec = AMP::LinearAlgebra::createVector(tempDofMap, temperatureVariable , true);
  AMP::LinearAlgebra::Vector::shared_ptr finalTempVec = initTempVec->cloneVector();

  initTempVec->setRandomValues();
  initTempVec->abs ( initTempVec );
  double initTempConst = input_db->getDoubleWithDefault("INIT_TEMP_CONST", 10.0);
  initTempVec->scale(initTempConst);
  initTempVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);

  bool setFinalTempEqualsInitialTemp = input_db->getBoolWithDefault("SET_FINAL_TEMP_EQUALS_INIT_TEMP", false);

  if(setFinalTempEqualsInitialTemp) {
    finalTempVec->copyVector(initTempVec);
  } else {
    finalTempVec->setRandomValues();
    double finalTempConst = input_db->getDoubleWithDefault("FINAL_TEMP_CONST", 12.0);
    finalTempVec->scale(finalTempConst);
  }
  finalTempVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);

  (boost::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(nonlinBvpOperator->
                                                                            getVolumeOperator()))->setReferenceTemperature(initTempVec);
  (boost::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(nonlinBvpOperator->
                                                                            getVolumeOperator()))->setVector(AMP::Operator::Mechanics::TEMPERATURE, finalTempVec); 

  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linBvpOperator =
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "linearMechanicsBVPOperator", input_db, elementPhysicsModel));

  //For RHS (Point Forces)
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletLoadVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "Load_Boundary", input_db, dummyModel));
  dirichletLoadVecOp->setVariable(displacementVariable);

  //For Initial-Guess
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletDispInVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "Displacement_Boundary",
          input_db,
          dummyModel));
  dirichletDispInVecOp->setVariable(displacementVariable);

  AMP::Discretization::DOFManager::shared_ptr dispDofMap = AMP::Discretization::simpleDOFManager::create(
      meshAdapter, AMP::Mesh::Vertex, 1, 3, true); 

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr mechNlSolVec = AMP::LinearAlgebra::createVector(dispDofMap, displacementVariable, true);
  AMP::LinearAlgebra::Vector::shared_ptr mechNlRhsVec = mechNlSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr mechNlResVec = mechNlSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr mechNlScaledRhsVec = mechNlSolVec->cloneVector();


  //Initial guess for NL solver must satisfy the displacement boundary conditions
  mechNlSolVec->setToScalar(0.0);
  dirichletDispInVecOp->apply(nullVec, nullVec, mechNlSolVec, 1.0, 0.0);
  mechNlSolVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);

  nonlinBvpOperator->apply(nullVec, mechNlSolVec, mechNlResVec, 1.0, 0.0);
  linBvpOperator->reset(nonlinBvpOperator->getJacobianParameters(mechNlSolVec));

  mechNlRhsVec->setToScalar(0.0);
  dirichletLoadVecOp->apply(nullVec, nullVec, mechNlRhsVec, 1.0, 0.0);

  double initSolNorm = mechNlSolVec->L2Norm();

  std::cout<<"Initial Solution Norm: "<<initSolNorm<<std::endl;

  boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 

  boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

  // ---- first initialize the preconditioner
  boost::shared_ptr<AMP::Database> pcSolver_db = linearSolver_db->getDatabase("Preconditioner"); 
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> pcSolverParams(new AMP::Solver::SolverStrategyParameters(pcSolver_db));
  pcSolverParams->d_pOperator = linBvpOperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(new AMP::Solver::TrilinosMLSolver(pcSolverParams));

  //HACK to prevent a double delete on Petsc Vec
  boost::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver;

  // initialize the linear solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));

  linearSolverParams->d_pOperator = linBvpOperator;
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pPreconditioner = pcSolver;

  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));

  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new
      AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm = globalComm;
  nonlinearSolverParams->d_pOperator = nonlinBvpOperator;
  nonlinearSolverParams->d_pKrylovSolver = linearSolver;
  nonlinearSolverParams->d_pInitialGuess = mechNlSolVec;

  nonlinearSolver.reset(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

  nonlinearSolver->setZeroInitialGuess(false);

#ifdef USE_EXT_SILO
  siloWriter->registerVector(mechNlSolVec, meshAdapter, AMP::Mesh::Vertex, "Solution" );
#endif

  for (int step=0;step<NumberOfLoadingSteps; step++)
  {
    AMP::pout << "########################################" << std::endl;
    AMP::pout << "The current loading step is " << (step+1) << std::endl;

    double scaleValue  = ((double)step+1.0)/NumberOfLoadingSteps;
    mechNlScaledRhsVec->scale(scaleValue, mechNlRhsVec);
    mechNlScaledRhsVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
    AMP::pout << "L2 Norm at loading step " << (step+1) << " is " << mechNlScaledRhsVec->L2Norm() << std::endl;

    nonlinBvpOperator->apply(mechNlScaledRhsVec, mechNlSolVec, mechNlResVec, 1.0, -1.0);
    double initialResidualNorm  = mechNlResVec->L2Norm();
    AMP::pout<<"Initial Residual Norm for loading step "<<(step+1)<<" is "<<initialResidualNorm<<std::endl;

    if(initialResidualNorm < 1.0e-2) {
      ut->passes("Nonlinear solve for current loading step");
    } else {
      AMP::pout<<"Starting Nonlinear Solve..."<<std::endl;
      nonlinearSolver->solve(mechNlScaledRhsVec, mechNlSolVec);

      nonlinBvpOperator->apply(mechNlScaledRhsVec, mechNlSolVec, mechNlResVec, 1.0, -1.0);
      double finalResidualNorm  = mechNlResVec->L2Norm();
      AMP::pout<<"Final Residual Norm for loading step "<<(step+1)<<" is "<<finalResidualNorm<<std::endl;

      if( finalResidualNorm > (1.0e-8*initialResidualNorm) ) {
        ut->failure("Nonlinear solve for current loading step");
      } else {
        ut->passes("Nonlinear solve for current loading step");
      }
    }

    boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
    boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperatorParameters> tmpParams(new
        AMP::Operator::MechanicsNonlinearFEOperatorParameters(tmp_db));
    (nonlinBvpOperator->getVolumeOperator())->reset(tmpParams);
    nonlinearSolver->setZeroInitialGuess(false);

#ifdef USE_EXT_SILO
    siloWriter->writeFile("FrozenTemp_NonlinearMechExample", step);
#endif
  }

  double finalSolNorm = mechNlSolVec->L2Norm();
  AMP::pout<<"Final Solution Norm: "<<finalSolNorm<<std::endl;

  ut->passes(exeName);

}


int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  try {
    myTest(&ut);
  } catch (std::exception &err) {
    std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
    ut.failure("ERROR: While testing");
  } catch( ... ) {
    std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
    ut.failure("ERROR: While testing");
  }

  ut.report();

  int num_failed = ut.NumFailGlobal();
  AMP::AMPManager::shutdown();
  return num_failed;
}   


