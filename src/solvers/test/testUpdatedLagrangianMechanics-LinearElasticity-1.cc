
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/libmesh/libMesh.h"

#include "discretization/simpleDOF_Manager.h"
#include "discretization/DOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "utils/Writer.h"

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

#include "utils/ReadTestMesh.h"
#include "libmesh/mesh_communication.h"

#include <iostream>
#include <string>

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "log_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  //Read the input file
  AMP::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  std::string mesh_file = input_db->getString("mesh_file");
  const unsigned int mesh_dim = 3;
  AMP::shared_ptr< ::Mesh > mesh(new ::Mesh(mesh_dim));
  AMP::readTestMesh(mesh_file, mesh);
  MeshCommunication().broadcast(*(mesh.get()));
  mesh->prepare_for_use(false);

  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::shared_ptr (
      new AMP::Mesh::libMesh (mesh, "TestMesh") );

  AMP_INSIST(input_db->keyExists("NumberOfLoadingSteps"), "Key ''NumberOfLoadingSteps'' is missing!");
  int NumberOfLoadingSteps = input_db->getInteger("NumberOfLoadingSteps");

  AMP_INSIST(input_db->keyExists("OutputFileName"), "Key ''OutputFileName'' is missing!");
  std::string outFileName = input_db->getString("OutputFileName");

  FILE* fp; 
  fp = fopen(outFileName.c_str(), "w");
  fprintf(fp, "clc; \n clear; \n A = zeros(24, 24); \n \n");

  //Create a nonlinear BVP operator for mechanics
  AMP_INSIST( input_db->keyExists("NonlinearMechanicsOperator"), "key missing!" );
  AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearMechanicsBVPoperator = AMP::dynamic_pointer_cast<
    AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(
          meshAdapter, "NonlinearMechanicsOperator", input_db ));
  AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> nonlinearMechanicsVolumeOperator = 
    AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(nonlinearMechanicsBVPoperator->getVolumeOperator());
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> mechanicsMaterialModel = nonlinearMechanicsVolumeOperator->getMaterialModel();

  //Create a Linear BVP operator for mechanics
  AMP_INSIST( input_db->keyExists("LinearMechanicsOperator"), "key missing!" );
  //AMP::shared_ptr<AMP::Database> linearMechanicsDatabase = input_db->getDatabase("LinearMechanicsOperator");
  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearMechanicsBVPoperator = AMP::dynamic_pointer_cast<
    AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(
          meshAdapter, "LinearMechanicsOperator", input_db, mechanicsMaterialModel));

  //Create the variables
  AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> mechanicsNonlinearVolumeOperator = 
    AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
        nonlinearMechanicsBVPoperator->getVolumeOperator());
  AMP::LinearAlgebra::Variable::shared_ptr dispVar = mechanicsNonlinearVolumeOperator->getOutputVariable();

  //For RHS (Point Forces)
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletLoadVecOp =
    AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
        AMP::Operator::OperatorBuilder::createOperator(meshAdapter, "Load_Boundary", input_db, dummyModel));
  dirichletLoadVecOp->setVariable(dispVar);

  AMP::Discretization::DOFManager::shared_ptr dofMap = AMP::Discretization::simpleDOFManager::create(
      meshAdapter, AMP::Mesh::Vertex, 1, 3, true); 

  //Create the vectors
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector(dofMap, dispVar, true);
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = solVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr resVec = solVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr scaledRhsVec = solVec->cloneVector();

  //Initial guess
  solVec->zero();

  //RHS
  rhsVec->zero();
  dirichletLoadVecOp->apply(nullVec, rhsVec);
  nonlinearMechanicsBVPoperator->modifyRHSvector(rhsVec);

  //We need to reset the linear operator before the solve since TrilinosML does
  //the factorization of the matrix during construction and so the matrix must
  //be correct before constructing the TrilinosML object.
  nonlinearMechanicsBVPoperator->apply( solVec, resVec);
  linearMechanicsBVPoperator->reset(nonlinearMechanicsBVPoperator->getJacobianParameters(solVec));

  double epsilon = 1.0e-13*(((linearMechanicsBVPoperator->getMatrix())->extractDiagonal())->L1Norm());

  AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
  AMP::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

  // ---- first initialize the preconditioner
  AMP::shared_ptr<AMP::Database> pcSolver_db = linearSolver_db->getDatabase("Preconditioner"); 
  AMP::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> pcSolverParams(new AMP::Solver::TrilinosMLSolverParameters(pcSolver_db));
  pcSolverParams->d_pOperator = linearMechanicsBVPoperator;
  AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(new AMP::Solver::TrilinosMLSolver(pcSolverParams));

  //HACK to prevent a double delete on Petsc Vec
  AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver;

  // initialize the linear solver
  AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));
  linearSolverParams->d_pOperator = linearMechanicsBVPoperator;
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pPreconditioner = pcSolver;
  AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));

  // initialize the nonlinear solver
  AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new
      AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));
  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm = globalComm;
  nonlinearSolverParams->d_pOperator = nonlinearMechanicsBVPoperator;
  nonlinearSolverParams->d_pKrylovSolver = linearSolver;
  nonlinearSolverParams->d_pInitialGuess = solVec;
  nonlinearSolver.reset(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

  nonlinearSolver->setZeroInitialGuess(false);
  /*
     FILE *fout1;
     fout1 = fopen("Loading_Loop.txt","w");
     */
  //double delta_displacement_shear = 0.00006;
  //double delta_displacement_axial = 0.00006;

  for (int step=0;step<NumberOfLoadingSteps; step++)
  {
    AMP::pout << "########################################" << std::endl;
    AMP::pout << "The current loading step is " << (step+1) << std::endl;

    nonlinearMechanicsBVPoperator->modifyInitialSolutionVector(solVec);

    double scaleValue;
    /*    double No5 = NumberOfLoadingSteps / 5;
          double N2o5 = (2 * NumberOfLoadingSteps) / 5;
          double N3o5 = (3 * NumberOfLoadingSteps) / 5;
          double N4o5 = (4 * NumberOfLoadingSteps) / 5;
          double N = NumberOfLoadingSteps;
          if(step < No5) {
          scaleValue = ((double)step+1.0) / ((double)No5);
          }
          if((step >= No5) && (step < N2o5)) {
          scaleValue = 1.0 - (((double)step + 1.0 - (double)No5) / ((double)No5));
          }
          if((step >= N2o5) && (step < N3o5)) {
          scaleValue = - (((double)step + 1.0 - (double)N2o5) / ((double)No5));
          }
          if((step >= N3o5) && (step < N4o5)) {
          scaleValue = -1.0 + (((double)step + 1.0 - (double)N3o5) / ((double)No5));
          }
          if((step >= N4o5) && (step < N)) {
          scaleValue = (((double)step + 1.0 - (double)N4o5) / ((double)No5));
          }
          fprintf(fout1,"%lf %lf\n",((double)step + 1.0),scaleValue);
          */
    /*    double No2 = NumberOfLoadingSteps / 2;
          double N = NumberOfLoadingSteps;
          if(step < No2) {
          scaleValue = ((double)step+1.0) / ((double)No2);
          }
          if((step >= No2) && (step < N)) {
          scaleValue = 1.0 - (((double)step + 1.0 - (double)No2) / ((double)No2));
          }
          fprintf(fout1,"%lf %lf\n",((double)step + 1.0),scaleValue);
          */
    //double No4 = NumberOfLoadingSteps / 4;
    //double No2 = NumberOfLoadingSteps / 2;
    //double N3o4 = (3 * NumberOfLoadingSteps) / 4;
    //double N = NumberOfLoadingSteps;

    scaleValue  = ((double)step+1.0)/NumberOfLoadingSteps;
    scaledRhsVec->scale(scaleValue, rhsVec);
    AMP::pout << "L2 Norm of RHS at loading step " << (step+1) << " is " << scaledRhsVec->L2Norm() << std::endl;

    nonlinearMechanicsBVPoperator->residual(scaledRhsVec, solVec, resVec);
    double initialResidualNorm  = resVec->L2Norm();
    AMP::pout<<"Initial Residual Norm for loading step "<<(step+1)<<" is "<<initialResidualNorm<<std::endl;

    nonlinearSolver->solve(scaledRhsVec, solVec);

    nonlinearMechanicsBVPoperator->residual(scaledRhsVec, solVec, resVec);
    double finalResidualNorm  = resVec->L2Norm();
    AMP::pout<<"Final Residual Norm for loading step "<<(step+1)<<" is "<<finalResidualNorm<<std::endl;

    if( finalResidualNorm > (1.0e-10*initialResidualNorm) ) {
      ut->failure("Nonlinear solve for current loading step");
    } else {
      ut->passes("Nonlinear solve for current loading step");
    }

    double finalSolNorm = solVec->L2Norm();

    AMP::pout<<"Final Solution Norm: "<<finalSolNorm<<std::endl;

    AMP::LinearAlgebra::Vector::shared_ptr mechUvec = solVec->select( AMP::LinearAlgebra::VS_Stride(0,3), "U" );
    AMP::LinearAlgebra::Vector::shared_ptr mechVvec = solVec->select( AMP::LinearAlgebra::VS_Stride(1,3), "V" );
    AMP::LinearAlgebra::Vector::shared_ptr mechWvec = solVec->select( AMP::LinearAlgebra::VS_Stride(2,3), "W" );

    double finalMaxU = mechUvec->maxNorm();
    double finalMaxV = mechVvec->maxNorm();
    double finalMaxW = mechWvec->maxNorm();

    AMP::pout<<"Maximum U displacement: "<<finalMaxU<<std::endl;
    AMP::pout<<"Maximum V displacement: "<<finalMaxV<<std::endl;
    AMP::pout<<"Maximum W displacement: "<<finalMaxW<<std::endl;

    AMP::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
    AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperatorParameters> tmpParams(new
        AMP::Operator::MechanicsNonlinearFEOperatorParameters(tmp_db));
    (nonlinearMechanicsBVPoperator->getVolumeOperator())->reset(tmpParams);

    /*    if(step < (No4 - 1)) {
          dirichletVectorCorrectionDatabase->putDouble("value_1_0", (((double)(step + 2))*delta_displacement_axial));
          dirichletVectorCorrectionDatabase->putDouble("value_2_0", 0.0);
          AMP::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> bndParams(new 
          AMP::Operator::DirichletVectorCorrectionParameters(dirichletVectorCorrectionDatabase));
          (nonlinearMechanicsBVPoperator->getBoundaryOperator())->reset(bndParams);
          fprintf(fout1,"%d %le %le\n",step,(((double)(step + 2))*delta_displacement_axial),0.0);
          }
          if((step >= (No4 - 1)) && (step < (No2 - 1))) {
          dirichletVectorCorrectionDatabase->putDouble("value_1_0", 0.6);
          dirichletVectorCorrectionDatabase->putDouble("value_2_0", (((double)(step + 2 - No4))*delta_displacement_shear));
          AMP::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> bndParams(new 
          AMP::Operator::DirichletVectorCorrectionParameters(dirichletVectorCorrectionDatabase));
          (nonlinearMechanicsBVPoperator->getBoundaryOperator())->reset(bndParams);
          fprintf(fout1,"%d %le %le\n",step,0.3,(((double)(step + 2 - No4))*delta_displacement_shear));
          }
          if((step >= (No2 - 1)) && (step < (N3o4 - 1))) {
          dirichletVectorCorrectionDatabase->putDouble("value_1_0", (0.6 - (((double)(step + 2 - No2))*delta_displacement_axial)));
          dirichletVectorCorrectionDatabase->putDouble("value_2_0", 0.6);
          AMP::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> bndParams(new 
          AMP::Operator::DirichletVectorCorrectionParameters(dirichletVectorCorrectionDatabase));
          (nonlinearMechanicsBVPoperator->getBoundaryOperator())->reset(bndParams);
          fprintf(fout1,"%d %le %le\n",step,(0.3 - (((double)(step + 2 - No2))*delta_displacement_axial)),0.3);
          }
          if((step >= (N3o4 - 1)) && (step < N)) {
          dirichletVectorCorrectionDatabase->putDouble("value_1_0", 0.0);
          dirichletVectorCorrectionDatabase->putDouble("value_2_0", (0.6 - (((double)(step + 2 - N3o4))*delta_displacement_shear)));
          AMP::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> bndParams(new 
          AMP::Operator::DirichletVectorCorrectionParameters(dirichletVectorCorrectionDatabase));
          (nonlinearMechanicsBVPoperator->getBoundaryOperator())->reset(bndParams);
          fprintf(fout1,"%d %le %le\n",step,0.0,(0.3 - (((double)(step + 2 - N3o4))*delta_displacement_shear)));
          }
          */    

    /*    if(step == 0) {
          dirichletVectorCorrectionDatabase->putDouble("value_1_0", 0.3);
          dirichletVectorCorrectionDatabase->putDouble("value_2_0", 0.3);
          AMP::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> bndParams(new 
          AMP::Operator::DirichletVectorCorrectionParameters(dirichletVectorCorrectionDatabase));
          (nonlinearMechanicsBVPoperator->getBoundaryOperator())->reset(bndParams);
          fprintf(fout1,"%d %le %le\n",step,(((double)(step + 2))*0.003),0.0);
          }
          if(step == 1) {
          dirichletVectorCorrectionDatabase->putDouble("value_1_0", 0.0);
          dirichletVectorCorrectionDatabase->putDouble("value_2_0", 0.3);
          AMP::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> bndParams(new 
          AMP::Operator::DirichletVectorCorrectionParameters(dirichletVectorCorrectionDatabase));
          (nonlinearMechanicsBVPoperator->getBoundaryOperator())->reset(bndParams);
          fprintf(fout1,"%d %le %le\n",step,0.3,(((double)(step + 2 - No4))*0.003));
          }
          if((step == 2) || (step == 3)) {
          dirichletVectorCorrectionDatabase->putDouble("value_1_0", 0.0);
          dirichletVectorCorrectionDatabase->putDouble("value_2_0", 0.0);
          AMP::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> bndParams(new 
          AMP::Operator::DirichletVectorCorrectionParameters(dirichletVectorCorrectionDatabase));
          (nonlinearMechanicsBVPoperator->getBoundaryOperator())->reset(bndParams);
          fprintf(fout1,"%d %le %le\n",step,(0.3 - (((double)(step + 2 - No2))*0.003)),0.3);
          }
          */    
    nonlinearSolver->setZeroInitialGuess(false);

    //std::cout<<solVec<<std::endl;

    AMP::shared_ptr<AMP::LinearAlgebra::Matrix> mechMat = linearMechanicsBVPoperator->getMatrix();

    for(int i = 0; i < 24; i++) {
      std::vector<unsigned int> matCols;
      std::vector<double> matVals;
      mechMat->getRowByGlobalID(i, matCols, matVals);
      for(unsigned int j = 0; j < matCols.size(); j++) {
        fprintf(fp, "A(%d, %d) = %.15f ; \n", (i + 1), (int)(matCols[j] + 1), matVals[j]);
      }//end for j
      fprintf(fp, "\n");
    }//end for i

    /*char num1[256];
      sprintf(num1,"%d",step);
      std::string number1 = num1;
      std::string fname = exeName + "_Stress_Strain_" + number1 + ".txt";

      AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(nonlinearMechanicsBVPoperator->getVolumeOperator())->printStressAndStrain(solVec, fname);
      */
  }

  //fclose(fout1);

  //AMP::pout<<solVec<<std::endl;

  AMP::pout<<"epsilon = "<<epsilon<<std::endl;

  //mechanicsNonlinearVolumeOperator->printStressAndStrain(solVec, output_file);

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit(new AMP::Mesh::initializeLibMesh(AMP_COMM_WORLD));

  AMP::UnitTest ut;

  std::vector<std::string> exeNames;
  //exeNames.push_back("testUpdatedLagrangianMechanics-LinearElasticity-1");
  //exeNames.push_back("testUpdatedLagrangianMechanics-LinearElasticity-1a");
  exeNames.push_back("testUpdatedLagrangianMechanics-LinearElasticity-1b");
  //exeNames.push_back("testUpdatedLagrangianMechanics-LinearElasticity-1c");

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

  libmeshInit.reset();
  AMP::AMPManager::shutdown();
  return num_failed;
}   



