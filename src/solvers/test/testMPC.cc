
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "ampmesh/MeshVariable.h"

#include "materials/Material.h"

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/ReadTestMesh.h"
#include "utils/WriteSolutionToFile.h"

#include "operators/PetscMatrixShellOperator.h"
#include "operators/ColumnOperator.h"
#include "operators/ContactResidualCorrection.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/mechanics/MechanicsLinearFEOperator.h"

#include "solvers/ColumnSolver.h"
#include "solvers/PetscKrylovSolverParameters.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/TrilinosMLSolver.h"
#include "solvers/ContactSolver.h"

#include "mesh_communication.h"

#include "MPCtestUtils.h"

extern "C"{
#include "petsc.h"
}

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);
  std::string mesh_file = input_db->getString("mesh_file");

  AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );

  const unsigned int mesh_dim = 3;
  boost::shared_ptr< ::Mesh > meshMaster(new ::Mesh(mesh_dim));
  boost::shared_ptr< ::Mesh > meshSlave(new ::Mesh(mesh_dim));

  bool binaryMeshes = input_db->getBool("BinaryMeshes");

  if(binaryMeshes) {
    AMP::readBinaryTestMesh(mesh_file, meshMaster);
    AMP::readBinaryTestMesh(mesh_file, meshSlave);
  } else {
    AMP::readTestMesh(mesh_file, meshMaster);
    AMP::readTestMesh(mesh_file, meshSlave);
  }

  MeshCommunication().broadcast(*(meshMaster.get()));
  MeshCommunication().broadcast(*(meshSlave.get()));

  meshMaster->prepare_for_use(false);
  meshSlave->prepare_for_use(false);

  AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter ( new AMP::Mesh::MeshManager::Adapter (meshMaster) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter ( new AMP::Mesh::MeshManager::Adapter (meshSlave) );

  manager->addMesh(masterMeshAdapter, "master");
  manager->addMesh(slaveMeshAdapter, "slave");

  double xOffset = input_db->getDouble("xOffset");
  double yOffset = input_db->getDouble("yOffset");
  double zOffset = input_db->getDouble("zOffset");
  slaveMeshAdapter->translate(xOffset, yOffset, zOffset);

  std::vector<unsigned int> masterNodes;
  std::vector<unsigned int> slaveNodes;

  getMasterAndSlaveNodes(masterNodes, slaveNodes, masterMeshAdapter, slaveMeshAdapter);

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> masterElementPhysicsModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> masterOperator = boost::dynamic_pointer_cast<
  AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(masterMeshAdapter,
										   "MasterOperator",
										   input_db,
										   masterElementPhysicsModel));

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel1;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> slaveOperator1 = boost::dynamic_pointer_cast<
    AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter,
										     "SlaveOperator",
										     input_db,
										     slaveElementPhysicsModel1));

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel2;
  boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> slaveOperator2 = boost::dynamic_pointer_cast<
    AMP::Operator::MechanicsLinearFEOperator>(AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter,
											     "VolumeOperator",
											     input_db,
											     slaveElementPhysicsModel2));

  AMP::LinearAlgebra::Variable::shared_ptr masterVar = masterOperator->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr slaveVar = slaveOperator1->getOutputVariable();

  boost::shared_ptr<AMP::Database> contactOpDatabase = input_db->getDatabase("ContactOperator");
  boost::shared_ptr<AMP::Operator::ContactResidualCorrectionParameters> contactParams(new 
      AMP::Operator::ContactResidualCorrectionParameters(contactOpDatabase));
  boost::shared_ptr<AMP::Operator::ContactResidualCorrection> contactOperator(new 
      AMP::Operator::ContactResidualCorrection(contactParams));   
  contactOperator->setMasterMesh(masterMeshAdapter);
  contactOperator->setSlaveMesh(slaveMeshAdapter);
  contactOperator->setMasterVariable(masterVar);
  contactOperator->setSlaveVariable(slaveVar);

  std::vector<std::vector<unsigned int> > dofs;
  for(size_t i = 0; i < masterNodes.size(); i++) {
    std::vector<unsigned int> tmpVec(3);
    tmpVec[0] = 0;
    tmpVec[1] = 1;
    tmpVec[2] = 2;
    dofs.push_back(tmpVec);
  }//end for i

  contactOperator->setMasterNodes(masterNodes);
  contactOperator->setSlaveNodes(slaveNodes);
  contactOperator->setDofs(dofs);

  boost::shared_ptr<AMP::Operator::OperatorParameters> dummyParams;
  boost::shared_ptr<AMP::Operator::ColumnOperator> columnOperator1(new AMP::Operator::ColumnOperator(dummyParams));
  columnOperator1->append(masterOperator);
  columnOperator1->append(slaveOperator1);
  columnOperator1->append(contactOperator);

  boost::shared_ptr<AMP::Operator::ColumnOperator> columnOperator2(new AMP::Operator::ColumnOperator(dummyParams));
  columnOperator2->append(masterOperator);
  columnOperator2->append(slaveOperator2);
  columnOperator2->append(contactOperator);

  AMP::LinearAlgebra::Variable::shared_ptr columnVar = contactOperator->getOutputVariable();

  boost::shared_ptr<AMP::Database> matrixShellDatabase = input_db->getDatabase("MatrixShellOperator");
  boost::shared_ptr<AMP::Operator::OperatorParameters> matrixShellParams(new
      AMP::Operator::OperatorParameters(matrixShellDatabase));
  boost::shared_ptr<AMP::Operator::PetscMatrixShellOperator> matrixShellOperator(new 
      AMP::Operator::PetscMatrixShellOperator(matrixShellParams));   

  int numMasterLocalNodes = masterMeshAdapter->numLocalNodes();
  int numSlaveLocalNodes = slaveMeshAdapter->numLocalNodes();
  int matLocalSize = (3*(numMasterLocalNodes + numSlaveLocalNodes));
  std::cout<<"numMasterNodes = "<<numMasterLocalNodes<<std::endl;
  std::cout<<"numSlaveNodes = "<<numSlaveLocalNodes<<std::endl;
  std::cout<<"matLocalSize = "<<matLocalSize<<std::endl;

  matrixShellOperator->setComm(globalComm);
  matrixShellOperator->setMatLocalRowSize(matLocalSize);
  matrixShellOperator->setMatLocalColumnSize(matLocalSize);
  matrixShellOperator->setOperator(columnOperator2); 

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOperator =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
									  AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter, "Load_Boundary", input_db, dummyModel));
  loadOperator->setVariable(slaveVar);

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr columnSolVec = manager->createVector(columnVar);
  AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec = manager->createVector(columnVar);

  columnSolVec->zero();
  columnRhsVec->zero();
  loadOperator->apply(nullVec, nullVec, columnRhsVec, 1.0, 0.0);

  boost::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase("LinearSolver"); 

  boost::shared_ptr<AMP::Database> columnPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::ColumnSolverParameters> columnPreconditionerParams(new
      AMP::Solver::ColumnSolverParameters(columnPreconditioner_db));
  columnPreconditionerParams->d_pOperator = columnOperator1;
  boost::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(new AMP::Solver::ColumnSolver(columnPreconditionerParams));

  bool useTrilinosPC = columnPreconditioner_db->getBool("useTrilinosPC");

  if(useTrilinosPC) {
    boost::shared_ptr<AMP::Database> masterPreconditioner_db = columnPreconditioner_db->getDatabase("MasterPreconditioner"); 
    boost::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> masterPreconditionerParams(new 
        AMP::Solver::TrilinosMLSolverParameters(masterPreconditioner_db));
    masterPreconditionerParams->d_pOperator = masterOperator;
    boost::shared_ptr<AMP::Solver::TrilinosMLSolver> masterPreconditioner(new AMP::Solver::TrilinosMLSolver(masterPreconditionerParams));
    columnPreconditioner->append(masterPreconditioner);

    boost::shared_ptr<AMP::Database> slavePreconditioner_db = columnPreconditioner_db->getDatabase("SlavePreconditioner"); 
    boost::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> slavePreconditionerParams(new 
        AMP::Solver::TrilinosMLSolverParameters(slavePreconditioner_db));
    slavePreconditionerParams->d_pOperator = slaveOperator1;
    boost::shared_ptr<AMP::Solver::TrilinosMLSolver> slavePreconditioner(new AMP::Solver::TrilinosMLSolver(slavePreconditionerParams));
    columnPreconditioner->append(slavePreconditioner);
  } else {
    boost::shared_ptr<AMP::Database> masterSolver_db = columnPreconditioner_db->getDatabase("MasterSolver"); 
    boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> masterSolverParams(new
        AMP::Solver::PetscKrylovSolverParameters(masterSolver_db));
    masterSolverParams->d_pOperator = masterOperator;
    masterSolverParams->d_comm = globalComm;
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> masterSolver(new AMP::Solver::PetscKrylovSolver(masterSolverParams));
    columnPreconditioner->append(masterSolver);

    boost::shared_ptr<AMP::Database> slaveSolver_db = columnPreconditioner_db->getDatabase("SlaveSolver"); 
    boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> slaveSolverParams(new
        AMP::Solver::PetscKrylovSolverParameters(slaveSolver_db));
    slaveSolverParams->d_pOperator = slaveOperator1;
    slaveSolverParams->d_comm = globalComm;
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> slaveSolver(new AMP::Solver::PetscKrylovSolver(slaveSolverParams));
    columnPreconditioner->append(slaveSolver);
  }

  boost::shared_ptr<AMP::Database> contactPreconditioner_db = columnPreconditioner_db->getDatabase("ContactPreconditioner"); 
  boost::shared_ptr<AMP::Solver::ContactSolverParameters> contactPreconditionerParams(new 
      AMP::Solver::ContactSolverParameters(contactPreconditioner_db));
  contactPreconditionerParams->d_pOperator = contactOperator;
  boost::shared_ptr<AMP::Solver::ContactSolver> contactPreconditioner(new AMP::Solver::ContactSolver(contactPreconditionerParams));
  columnPreconditioner->append(contactPreconditioner);

  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));
  linearSolverParams->d_pOperator = matrixShellOperator;
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pPreconditioner = columnPreconditioner;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));
  linearSolver->setZeroInitialGuess(true);

  linearSolver->solve(columnRhsVec, columnSolVec);

  AMP::LinearAlgebra::Vector::shared_ptr masterSolVec = columnSolVec->subsetVectorForVariable(masterVar);
  AMP::LinearAlgebra::Vector::shared_ptr slaveSolVec = columnSolVec->subsetVectorForVariable(slaveVar);

  printSolution(masterMeshAdapter, masterSolVec, (exeName + "-master"));
  printSolution(slaveMeshAdapter, slaveSolVec, (exeName + "-slave"));

  ut->passes(exeName);
}


int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::vector<std::string> exeNames;

  if(argc == 1) {
    exeNames.push_back("testMPC-1");
    exeNames.push_back("testMPC-2");
    exeNames.push_back("testMPC-3");
    exeNames.push_back("testMPC-4");
    exeNames.push_back("testMPC-5");
  } else {
    for(int i = 1; i < argc; i++) {
      char inpName[100];
      sprintf(inpName, "testMPC-%s", argv[i]);
      exeNames.push_back(inpName);
    }//end for i
  }

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

