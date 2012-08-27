#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/SiloIO.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"

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

  // Get the Mesh database and create the mesh parameters
  boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
  boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
  params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));

  // Create the meshes from the input database
  AMP::Mesh::Mesh::shared_ptr  mesh = AMP::Mesh::Mesh::buildMesh(params);

  // Create the DOFManagers
  AMP::Discretization::DOFManager::shared_ptr NodalVectorDOF = 
    AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,3);

  AMP_INSIST(input_db->keyExists("NumberOfLoadingSteps"), "Key ''NumberOfLoadingSteps'' is missing!");
  int NumberOfLoadingSteps = input_db->getInteger("NumberOfLoadingSteps");

  bool ExtractData = input_db->getBoolWithDefault("ExtractStressStrainData", false);
  FILE *fout123;
  std::string ss_file = exeName + "_UniaxialTmperatureDisplacement.txt";
  fout123 = fopen(ss_file.c_str(),"w");

  //Create a nonlinear BVP operator for mechanics
  AMP_INSIST( input_db->keyExists("NonlinearMechanicsOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> mechanicsMaterialModel;
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearMechanicsBVPoperator = boost::dynamic_pointer_cast<
  AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(mesh,
										      "NonlinearMechanicsOperator",
										      input_db,
										      mechanicsMaterialModel));

  //Create a Linear BVP operator for mechanics
  AMP_INSIST( input_db->keyExists("LinearMechanicsOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearMechanicsBVPoperator = boost::dynamic_pointer_cast<
  AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(mesh,
										   "LinearMechanicsOperator",
										   input_db,
										   mechanicsMaterialModel));

  //Create the variables
  boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> mechanicsNonlinearVolumeOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
        nonlinearMechanicsBVPoperator->getVolumeOperator());

  boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> multivariable = boost::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(
      mechanicsNonlinearVolumeOperator->getInputVariable()); 
  AMP::LinearAlgebra::Variable::shared_ptr dispVar = multivariable->getVariable(AMP::Operator::Mechanics::DISPLACEMENT); 
  AMP::LinearAlgebra::Variable::shared_ptr tempVar =  multivariable->getVariable(AMP::Operator::Mechanics::TEMPERATURE); 
  AMP::LinearAlgebra::Variable::shared_ptr burnVar =  multivariable->getVariable(AMP::Operator::Mechanics::BURNUP); 

  //boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> mechanicsLinearVolumeOperator = 
  //  boost::dynamic_pointer_cast<AMP::Operator::MechanicsLinearFEOperator>(
  //      linearMechanicsBVPoperator->getVolumeOperator());

  //For RHS (Point Forces)
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletLoadVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(mesh,
															 "Load_Boundary",
															 input_db,
															 dummyModel));
  dirichletLoadVecOp->setVariable(dispVar);

  //Create the vectors
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( NodalVectorDOF, dispVar );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = solVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr resVec = solVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr tempVec = AMP::LinearAlgebra::createVector( NodalVectorDOF, tempVar );
  AMP::LinearAlgebra::Vector::shared_ptr tempVecRef = tempVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr burnVec = AMP::LinearAlgebra::createVector( NodalVectorDOF, burnVar );

  //Initial guess
  solVec->zero();
  nonlinearMechanicsBVPoperator->modifyInitialSolutionVector(solVec);

  //RHS
  rhsVec->zero();
  dirichletLoadVecOp->apply(nullVec, nullVec, rhsVec, 1.0, 0.0);
  nonlinearMechanicsBVPoperator->modifyRHSvector(rhsVec);

  // Create the silo writer and register the data
  #ifdef USES_SILO
    AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO);
    siloWriter->registerVector( solVec, mesh, AMP::Mesh::Vertex, "Solution_Vector" );
  #endif

  //Adding the Temperature and Burnup
  tempVecRef->setToScalar(301.0);
  tempVec->setToScalar(301.0);
  burnVec->setToScalar(10.0);

  mechanicsNonlinearVolumeOperator->setReferenceTemperature(tempVecRef);
  mechanicsNonlinearVolumeOperator->setVector(AMP::Operator::Mechanics::TEMPERATURE, tempVec);
  mechanicsNonlinearVolumeOperator->setVector(AMP::Operator::Mechanics::BURNUP, burnVec);

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

  //double TotalLoadingSteps = NumberOfLoadingSteps / 4;
  //double TotalUnloadingSteps = NumberOfLoadingSteps - TotalLoadingSteps;

  if(ExtractData) {
    fprintf(fout123,"%lf %lf %lf %lf\n",301.0,0.0,0.0,0.0);
  }

  for (int step=0;step<NumberOfLoadingSteps; step++)
  {
    AMP::pout << "########################################" << std::endl;
    AMP::pout << "The current loading step is " << (step+1) << std::endl;

    double finalTemperature = 301.0 + (((double)(step + 1)) * 200.0);
    tempVec->setToScalar(finalTemperature);

    nonlinearMechanicsBVPoperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);
    double initialResidualNorm  = resVec->L2Norm();
    AMP::pout<<"Initial Residual Norm for loading step "<<(step+1)<<" is "<<initialResidualNorm<<std::endl;

    nonlinearSolver->solve(rhsVec, solVec);

    nonlinearMechanicsBVPoperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);
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

    //double prev_stress, prev_strain, slope;
    if(ExtractData) {
      fprintf(fout123,"%lf %lf %lf %lf\n",finalTemperature,finalMaxU,finalMaxV,finalMaxW);
    }
  }

  AMP::pout<<"epsilon = "<<epsilon<<std::endl;

  mechanicsNonlinearVolumeOperator->printStressAndStrain(solVec, output_file);

#ifdef USES_SILO
    
    siloWriter->writeFile(exeName, 1);
#endif

  ut->passes(exeName);
  fclose(fout123);
}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    //exeNames.push_back("testFixedBeam-ThermalExpansion-1");
    //exeNames.push_back("testFixedBeam-ThermalExpansion-2");
    exeNames.push_back("testFixedBeam-ThermalExpansion-3");

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



