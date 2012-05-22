#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/SiloIO.h"
#include "ampmesh/libmesh/libMesh.h"
#include "vectors/VectorBuilder.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "operators/mechanics/ThermalStrainMaterialModel.h"
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

#include "solvers/PetscKrylovSolverParameters.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/PetscSNESSolverParameters.h"
#include "solvers/PetscSNESSolver.h"
#include "solvers/TrilinosMLSolver.h"

#include "ReadTestMesh.h"
#include "mesh_communication.h"


#include <iostream>
#include <string>
#include <math.h>

#include "boost/shared_ptr.hpp"



void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string output_file = "output_" + exeName + ".txt";
  std::string log_file = "log_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  //Read the input file
  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

//--------------------------------------------------
//   Create the Mesh.
//--------------------------------------------------
  boost::shared_ptr<AMP::Mesh::initializeLibMesh>  libmeshInit(new AMP::Mesh::initializeLibMesh(AMP::AMP_MPI(AMP_COMM_WORLD)));

  std::string mesh_file = input_db->getString("mesh_file");
  const unsigned int mesh_dim = 3;
  boost::shared_ptr< ::Mesh > mesh(new ::Mesh(mesh_dim));
  AMP::readTestMesh(mesh_file, mesh);
  MeshCommunication().broadcast(*(mesh.get()));
  mesh->prepare_for_use(false);
  AMP::Mesh::Mesh::shared_ptr meshAdapter ( new AMP::Mesh::libMesh(mesh,"cook") );
//--------------------------------------------------

  AMP_INSIST(input_db->keyExists("NumberOfLoadingSteps"), "Key ''NumberOfLoadingSteps'' is missing!");
  int NumberOfLoadingSteps = input_db->getInteger("NumberOfLoadingSteps");

  bool ExtractData = input_db->getBoolWithDefault("ExtractStressStrainData", false);
  FILE *fout123;
  std::string ss_file = exeName + "_UniaxialStressStrain.txt";
  fout123 = fopen(ss_file.c_str(),"w");

  //Create a nonlinear BVP operator for mechanics
  AMP_INSIST( input_db->keyExists("NonlinearMechanicsOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> mechanicsMaterialModel;
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearMechanicsBVPoperator = boost::dynamic_pointer_cast<
  AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
										      "NonlinearMechanicsOperator",
										      input_db,
										      mechanicsMaterialModel));

  //Create a Linear BVP operator for mechanics
  AMP_INSIST( input_db->keyExists("LinearMechanicsOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearMechanicsBVPoperator = boost::dynamic_pointer_cast<
    AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
										     "LinearMechanicsOperator",
										     input_db,
										     mechanicsMaterialModel));

  //Create the variables
  boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> mechanicsNonlinearVolumeOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
        nonlinearMechanicsBVPoperator->getVolumeOperator());
  AMP::LinearAlgebra::Variable::shared_ptr dispVar = mechanicsNonlinearVolumeOperator->getOutputVariable();

  //boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> mechanicsLinearVolumeOperator = 
  //  boost::dynamic_pointer_cast<AMP::Operator::MechanicsLinearFEOperator>(
  //      linearMechanicsBVPoperator->getVolumeOperator());

  //For RHS (Point Forces)
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletLoadVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
															 "Load_Boundary",
															 input_db,
															 dummyModel));
  dirichletLoadVecOp->setVariable(dispVar);

//--------------------------------------------------
// Create a DOF manager for a nodal vector 
//--------------------------------------------------
  int DOFsPerNode = 3;
  int nodalGhostWidth = 1;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr nodalDofMap      = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth,      DOFsPerNode,    split);
//--------------------------------------------------

  //Create the vectors
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( nodalDofMap, dispVar );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = solVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr resVec = solVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr scaledRhsVec = solVec->cloneVector();

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

  int NumberOfLoops = 2;
  int TotalLoadingSteps = 4;
  NumberOfLoadingSteps = NumberOfLoops * (4 * TotalLoadingSteps);
  double AngleIncrement = (11.0 / 7.0) * (1.0 / ((double)TotalLoadingSteps));

  for (int step=0;step<NumberOfLoadingSteps; step++)
  {
    AMP::pout << "########################################" << std::endl;
    AMP::pout << "The current loading step is " << (step+1) << std::endl;

    double scaleValue;
    scaleValue = sin(((double)step) * AngleIncrement);
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

    char num1[256];
    sprintf(num1,"%d",step);
    std::string number1 = num1;
    std::string fname = exeName + "_Stress_Strain_" + number1 + ".txt";
    
    boost::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(nonlinearMechanicsBVPoperator->getVolumeOperator())->printStressAndStrain(solVec, fname);

    double prev_stress = 1.0, prev_strain = 1.0, slope = 1.0;
    if(ExtractData) {
      FILE *fin; 
      fin = fopen(fname.c_str(),"r");
      double coord[3], stress1[6], strain1[6];
      for(int ijk = 0; ijk < 8; ijk++) {
        for(int klm = 0; klm < 3; klm++) {
          fscanf(fin,"%lf",&coord[klm]);
        }
        for(int klm = 0; klm < 6; klm++) {
          fscanf(fin,"%lf",&stress1[klm]);
        }
        for(int klm = 0; klm < 6; klm++) {
          fscanf(fin,"%lf",&strain1[klm]);
        }
        if(ijk == 7) {
          if(step == 0) {
            slope = 0.0;
          } else {
            slope = (stress1[2] - prev_stress) / (strain1[2] - prev_strain);
          }
          fprintf(fout123,"%lf %lf %lf\n",strain1[2],stress1[2],slope);
          prev_stress = stress1[2];
          prev_strain = strain1[2];
        }
      }
      fclose(fin);
    }
  }

  AMP::pout<<"epsilon = "<<epsilon<<std::endl;

  mechanicsNonlinearVolumeOperator->printStressAndStrain(solVec, output_file);

  ut->passes(exeName);
  fclose(fout123);
}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

  std::vector<std::string> exeNames;
  //exeNames.push_back("testElementLevel-IsotropicElasticity");
  //exeNames.push_back("testElementLevel-VonMisesIsotropicHardeningPlasticity");
  exeNames.push_back("testElementLevel-VonMisesKinematicHardeningPlasticity");

    for(unsigned int i = 0; i < exeNames.size(); i++) {
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



