
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

#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"

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

void myTest(AMP::UnitTest *ut, std::string exeName)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

#ifdef USE_EXT_SILO
  // Create the silo writer and register the data
  AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
#endif

  AMP::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase("Mesh");
  AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);

  AMP_INSIST(input_db->keyExists("NumberOfLoadingSteps"), "Key ''NumberOfLoadingSteps'' is missing!");
  int NumberOfLoadingSteps = input_db->getInteger("NumberOfLoadingSteps");

  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinBvpOperator = 
    AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "nonlinearMechanicsBVPOperator", input_db, elementPhysicsModel));

  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linBvpOperator =
    AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "linearMechanicsBVPOperator", input_db, elementPhysicsModel));

  AMP::LinearAlgebra::Variable::shared_ptr displacementVariable = nonlinBvpOperator->getOutputVariable();

  //For RHS (Point Forces)
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletLoadVecOp =
    AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
        AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "Load_Boundary", input_db, dummyModel));
  dirichletLoadVecOp->setVariable(displacementVariable);

  //For Initial-Guess
  AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletDispInVecOp =
    AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
        AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "Displacement_Boundary", input_db, dummyModel));
  dirichletDispInVecOp->setVariable(displacementVariable);

  AMP::Discretization::DOFManager::shared_ptr dofMap = AMP::Discretization::simpleDOFManager::create(
      meshAdapter, AMP::Mesh::Vertex, 1, 3, true); 

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr mechNlSolVec = AMP::LinearAlgebra::createVector(dofMap, displacementVariable, true);
  AMP::LinearAlgebra::Vector::shared_ptr mechNlRhsVec = mechNlSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr mechNlResVec = mechNlSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr mechNlScaledRhsVec = mechNlSolVec->cloneVector();

  //Initial guess for NL solver must satisfy the displacement boundary conditions
  mechNlSolVec->setToScalar(0.0);
  dirichletDispInVecOp->apply( nullVec, mechNlSolVec);
  mechNlSolVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);

  nonlinBvpOperator->apply(mechNlSolVec, mechNlResVec);
  linBvpOperator->reset(nonlinBvpOperator->getJacobianParameters(mechNlSolVec));

  //Point forces
  mechNlRhsVec->setToScalar(0.0);
  dirichletLoadVecOp->apply(nullVec, mechNlRhsVec);

  AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
  AMP::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

  // ---- first initialize the preconditioner
  AMP::shared_ptr<AMP::Database> pcSolver_db = linearSolver_db->getDatabase("Preconditioner"); 
  AMP::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> pcSolverParams(new AMP::Solver::TrilinosMLSolverParameters(pcSolver_db));
  pcSolverParams->d_pOperator = linBvpOperator;
  AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(new AMP::Solver::TrilinosMLSolver(pcSolverParams));

  //HACK to prevent a double delete on Petsc Vec
  AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver;

  // initialize the linear solver
  AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));
  linearSolverParams->d_pOperator = linBvpOperator;
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pPreconditioner = pcSolver;
  AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));

  // initialize the nonlinear solver
  AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new
      AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));
  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm = globalComm;
  nonlinearSolverParams->d_pOperator = nonlinBvpOperator;
  nonlinearSolverParams->d_pKrylovSolver = linearSolver;
  nonlinearSolverParams->d_pInitialGuess = mechNlSolVec;
  nonlinearSolver.reset(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

  nonlinearSolver->setZeroInitialGuess(false);

  for (int step=0;step<NumberOfLoadingSteps; step++)
  {
    AMP::pout << "########################################" << std::endl;
    AMP::pout << "The current loading step is " << (step+1) << std::endl;

    double scaleValue  = ((double)step+1.0)/NumberOfLoadingSteps;
    mechNlScaledRhsVec->scale(scaleValue, mechNlRhsVec);
    mechNlScaledRhsVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
    AMP::pout << "L2 Norm of RHS at loading step " << (step+1) << " is " << mechNlScaledRhsVec->L2Norm() << std::endl;

    nonlinBvpOperator->residual(mechNlScaledRhsVec, mechNlSolVec, mechNlResVec);
    double initialResidualNorm  = mechNlResVec->L2Norm();
    AMP::pout<<"Initial Residual Norm for loading step "<<(step+1)<<" is "<<initialResidualNorm<<std::endl;

    AMP::pout<<"Starting Nonlinear Solve..."<<std::endl;
    nonlinearSolver->solve(mechNlScaledRhsVec, mechNlSolVec);

    nonlinBvpOperator->residual(mechNlScaledRhsVec, mechNlSolVec, mechNlResVec);
    double finalResidualNorm  = mechNlResVec->L2Norm();
    AMP::pout<<"Final Residual Norm for loading step "<<(step+1)<<" is "<<finalResidualNorm<<std::endl;

    if( finalResidualNorm > (1.0e-10*initialResidualNorm) ) {
      ut->failure("Nonlinear solve for current loading step");
    } else {
      ut->passes("Nonlinear solve for current loading step");
    }

    double finalSolNorm = mechNlSolVec->L2Norm();

    AMP::pout<<"Final Solution Norm: "<<finalSolNorm<<std::endl;

    AMP::LinearAlgebra::Vector::shared_ptr mechUvec = mechNlSolVec->select( AMP::LinearAlgebra::VS_Stride(0,3), "U" );
    AMP::LinearAlgebra::Vector::shared_ptr mechVvec = mechNlSolVec->select( AMP::LinearAlgebra::VS_Stride(1,3), "V" );
    AMP::LinearAlgebra::Vector::shared_ptr mechWvec = mechNlSolVec->select( AMP::LinearAlgebra::VS_Stride(2,3), "W" );

    double finalMaxU = mechUvec->maxNorm();
    double finalMaxV = mechVvec->maxNorm();
    double finalMaxW = mechWvec->maxNorm();

    AMP::pout<<"Maximum U displacement: "<<finalMaxU<<std::endl;
    AMP::pout<<"Maximum V displacement: "<<finalMaxV<<std::endl;
    AMP::pout<<"Maximum W displacement: "<<finalMaxW<<std::endl;

    AMP::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
    AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperatorParameters> tmpParams(new
        AMP::Operator::MechanicsNonlinearFEOperatorParameters(tmp_db));
    (nonlinBvpOperator->getVolumeOperator())->reset(tmpParams);
    nonlinearSolver->setZeroInitialGuess(false);

    meshAdapter->displaceMesh(mechNlSolVec);

#ifdef USE_EXT_SILO
    siloWriter->registerVector(mechNlSolVec, meshAdapter, AMP::Mesh::Vertex, "Solution" );
    char outFileName[256];
    sprintf(outFileName, "LoadPrescribed-DeformedPlateWithHole-LinearElasticity_%d", step);
    siloWriter->writeFile(outFileName, 0);
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

  std::vector<std::string> exeNames;
  exeNames.push_back("testPetscSNESSolver-NonlinearMechanics-PlateWithHole-1");
  exeNames.push_back("testPetscSNESSolver-LU-NonlinearMechanics-1-normal");
  exeNames.push_back("testPetscSNESSolver-ML-NonlinearMechanics-1-normal");
  exeNames.push_back("testPetscSNESSolver-LU-NonlinearMechanics-1-reduced");
  exeNames.push_back("testPetscSNESSolver-ML-NonlinearMechanics-1-reduced");

  for(size_t i = 0; i < exeNames.size(); i++) {
    try {
      myTest(&ut, exeNames[i]);
    } catch (std::exception &err) {
      std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
      ut.failure("ERROR: While testing");
    } catch( ... ) {
      std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
      ut.failure("ERROR: While testing");
    }
  }

  ut.report();

  int num_failed = ut.NumFailGlobal();
  AMP::AMPManager::shutdown();
  return num_failed;
}   


