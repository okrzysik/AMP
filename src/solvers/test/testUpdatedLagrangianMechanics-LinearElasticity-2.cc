#include "utils/AMPManager.h"
#include "utils/InputManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshVariable.h"
#include "ampmesh/SiloIO.h"

#include "operators/mechanics/ThermalStrainMaterialModel.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"

#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "solvers/PetscKrylovSolverParameters.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/PetscSNESSolverParameters.h"
#include "solvers/PetscSNESSolver.h"
#include "solvers/TrilinosMLSolver.h"

#include "ReadTestMesh.h"
#include "mesh_communication.h"


#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"



void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  //std::string output_file = "output_" + exeName + ".txt";
  std::string log_file = "log_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  //Read the input file
  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter;
  //meshAdapter = manager->getMesh ( "cylinder" );
  meshAdapter = manager->getMesh ( "brick" );
  
/*  std::string mesh_file = input_db->getString("mesh_file");
  const unsigned int mesh_dim = 3;
  boost::shared_ptr< ::Mesh > mesh(new ::Mesh(mesh_dim));
  AMP::readTestMesh(mesh_file, mesh);
  MeshCommunication().broadcast(*(mesh.get()));
  mesh->prepare_for_use(false);
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter ( new AMP::Mesh::MeshManager::Adapter (mesh) );
  manager->addMesh(meshAdapter, "cook");
*/
  AMP_INSIST(input_db->keyExists("NumberOfLoadingSteps"), "Key ''NumberOfLoadingSteps'' is missing!");
  int NumberOfLoadingSteps = input_db->getInteger("NumberOfLoadingSteps");

/*  AMP_INSIST(input_db->keyExists("OutputFileName"), "Key ''OutputFileName'' is missing!");
  std::string outFileName = input_db->getString("OutputFileName");

  FILE* fp; 
  fp = fopen(outFileName.c_str(), "w");
  fprintf(fp, "clc; \n clear; \n A = zeros(24, 24); \n \n");
*/
  //Create a nonlinear BVP operator for mechanics
  AMP_INSIST( input_db->keyExists("NonlinearMechanicsOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> mechanicsMaterialModel;
  //boost::shared_ptr<AMP::Database> nonlinearMechanicsDatabase = input_db->getDatabase("NonlinearMechanicsOperator");
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearMechanicsBVPoperator = boost::dynamic_pointer_cast<
  AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
										      "NonlinearMechanicsOperator",
										      input_db,
										      mechanicsMaterialModel));
  
  (boost::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(nonlinearMechanicsBVPoperator->getVolumeOperator()))->init();

  //Create a Linear BVP operator for mechanics
  AMP_INSIST( input_db->keyExists("LinearMechanicsOperator"), "key missing!" );
  //boost::shared_ptr<AMP::Database> linearMechanicsDatabase = input_db->getDatabase("LinearMechanicsOperator");
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearMechanicsBVPoperator = boost::dynamic_pointer_cast<
    AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
										     "LinearMechanicsOperator",
										     input_db,
										     mechanicsMaterialModel));

  //Create the variables
  boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> mechanicsNonlinearVolumeOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
        nonlinearMechanicsBVPoperator->getVolumeOperator());
  AMP::LinearAlgebra::Variable::shared_ptr dispVar = mechanicsNonlinearVolumeOperator->getInputVariable(AMP::Operator::Mechanics::DISPLACEMENT);

  //boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> mechanicsLinearVolumeOperator = 
  //  boost::dynamic_pointer_cast<AMP::Operator::MechanicsLinearFEOperator>(
  //      linearMechanicsBVPoperator->getVolumeOperator());

  //For RHS (Point Forces)
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  //boost::shared_ptr<AMP::Database> load_db = input_db->getDatabase("Load_Boundary");
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletLoadVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
															 "Load_Boundary",
															 input_db,
															 dummyModel));
  dirichletLoadVecOp->setVariable(dispVar);

  //Create the vectors
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr solVec = meshAdapter->createVector( dispVar );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = meshAdapter->createVector( dispVar );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = meshAdapter->createVector( dispVar );
  AMP::LinearAlgebra::Vector::shared_ptr scaledRhsVec = meshAdapter->createVector( dispVar );

  //Initial guess
  solVec->zero();
  nonlinearMechanicsBVPoperator->modifyInitialSolutionVector(solVec);

  //RHS
  rhsVec->zero();
  dirichletLoadVecOp->apply(nullVec, nullVec, rhsVec, 1.0, 0.0);
  nonlinearMechanicsBVPoperator->modifyRHSvector(rhsVec);

  //We need to reset the linear operator before the solve since TrilinosML does
  //the factorization of the matrix during construction and so the matrix must
  //be correct before constructing the TrilinosML object.
  nonlinearMechanicsBVPoperator->apply(nullVec, solVec, resVec, 1.0, 0.0);
  linearMechanicsBVPoperator->reset(nonlinearMechanicsBVPoperator->getJacobianParameters(solVec));

  double epsilon = 1.0e-13*(((linearMechanicsBVPoperator->getMatrix())->extractDiagonal())->L1Norm());

  boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
  boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

  // ---- first initialize the preconditioner
  boost::shared_ptr<AMP::Database> pcSolver_db = linearSolver_db->getDatabase("Preconditioner"); 
  boost::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> pcSolverParams(new AMP::Solver::TrilinosMLSolverParameters(pcSolver_db));
  pcSolverParams->d_pOperator = linearMechanicsBVPoperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(new AMP::Solver::TrilinosMLSolver(pcSolverParams));

  //HACK to prevent a double delete on Petsc Vec
  boost::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver;

  // initialize the linear solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));
  linearSolverParams->d_pOperator = linearMechanicsBVPoperator;
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pPreconditioner = pcSolver;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));

  // initialize the nonlinear solver
  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new
      AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));
  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm = globalComm;
  nonlinearSolverParams->d_pOperator = nonlinearMechanicsBVPoperator;
  nonlinearSolverParams->d_pKrylovSolver = linearSolver;
  nonlinearSolverParams->d_pInitialGuess = solVec;
  nonlinearSolver.reset(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

  nonlinearSolver->setZeroInitialGuess(false);

  int loadingSubSteps = 1;
  double scaleValue;

  for (int step=0;step<NumberOfLoadingSteps; step++)
  {
    AMP::pout << "########################################" << std::endl;
    AMP::pout << "The current loading step is " << (step+1) << std::endl;

    if(step > 3) {
      loadingSubSteps = 1;
    }
    for(int subStep = 0; subStep < loadingSubSteps; subStep++) {
      if(step <= 3) {
        scaleValue = ((double)step+1.0)/NumberOfLoadingSteps;
      } else {
        scaleValue = (((double)step) / NumberOfLoadingSteps) + (((double)subStep + 1.0) / (NumberOfLoadingSteps * loadingSubSteps));
      }

      AMP::pout << "########################################" << std::endl;
      AMP::pout << "The current loading sub step is " << (subStep+1) << " and scaleValue = " << scaleValue << std::endl;

      scaledRhsVec->scale(scaleValue, rhsVec);
      AMP::pout << "L2 Norm of RHS at loading step " << (step+1) << " is " << scaledRhsVec->L2Norm() << std::endl;

      nonlinearMechanicsBVPoperator->apply(scaledRhsVec, solVec, resVec, 1.0, -1.0);
      double initialResidualNorm  = resVec->L2Norm();
      AMP::pout<<"Initial Residual Norm for loading step "<<(step+1)<<" is "<<initialResidualNorm<<std::endl;

      nonlinearSolver->solve(scaledRhsVec, solVec);

      nonlinearMechanicsBVPoperator->apply(scaledRhsVec, solVec, resVec, 1.0, -1.0);
      double finalResidualNorm  = resVec->L2Norm();
      AMP::pout<<"Final Residual Norm for loading step "<<(step+1)<<" is "<<finalResidualNorm<<std::endl;

      if( finalResidualNorm > (1.0e-10*initialResidualNorm) ) {
        ut->failure("Nonlinear solve for current loading step");
      } else {
        ut->passes("Nonlinear solve for current loading step");
      }

      double finalSolNorm = solVec->L2Norm();

      AMP::pout<<"Final Solution Norm: "<<finalSolNorm<<std::endl;

      AMP::LinearAlgebra::Vector::shared_ptr mechUvec = solVec->select( AMP::LinearAlgebra::VS_Stride("U", 0, 3) , "U" );
      AMP::LinearAlgebra::Vector::shared_ptr mechVvec = solVec->select( AMP::LinearAlgebra::VS_Stride("V", 1, 3) , "V" );
      AMP::LinearAlgebra::Vector::shared_ptr mechWvec = solVec->select( AMP::LinearAlgebra::VS_Stride("W", 2, 3) , "W" );

      double finalMaxU = mechUvec->maxNorm();
      double finalMaxV = mechVvec->maxNorm();
      double finalMaxW = mechWvec->maxNorm();

      AMP::pout<<"Maximum U displacement: "<<finalMaxU<<std::endl;
      AMP::pout<<"Maximum V displacement: "<<finalMaxV<<std::endl;
      AMP::pout<<"Maximum W displacement: "<<finalMaxW<<std::endl;

      boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
      boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperatorParameters> tmpParams(new
          AMP::Operator::MechanicsNonlinearFEOperatorParameters(tmp_db));
      (nonlinearMechanicsBVPoperator->getVolumeOperator())->reset(tmpParams);
      nonlinearSolver->setZeroInitialGuess(false);
    }

    //std::cout<<solVec<<std::endl;

/*    boost::shared_ptr<AMP::LinearAlgebra::Matrix> mechMat = linearMechanicsBVPoperator->getMatrix();

    for(int i = 0; i < 24; i++) {
      std::vector<unsigned int> matCols;
      std::vector<double> matVals;
      mechMat->getRowByGlobalID(i, matCols, matVals);
      for(unsigned int j = 0; j < matCols.size(); j++) {
        fprintf(fp, "A(%d, %d) = %.15lf ; \n", (i + 1), (matCols[j] + 1), matVals[j]);
      }//end for j
      fprintf(fp, "\n");
    }//end for i
*/
/*    char num1[256];
    sprintf(num1,"%d",step);
    std::string number1 = num1;
    std::string fname = exeName + "_Stress_Strain_" + number1 + ".txt";
    
    boost::dynamic_pointer_cast<AMP::MechanicsNonlinearFEOperator>(nonlinearMechanicsBVPoperator->getVolumeOperator())->printStressAndStrain(solVec, fname);
*/    
  }

  AMP::pout<<"epsilon = "<<epsilon<<std::endl;

  //mechanicsNonlinearVolumeOperator->printStressAndStrain(solVec, output_file);

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    //exeNames.push_back("testUpdatedLagrangianMechanics-LinearElasticity-2");
    exeNames.push_back("testUpdatedLagrangianMechanics-LinearElasticity-3");

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



