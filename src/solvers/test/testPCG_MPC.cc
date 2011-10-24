
#include "utils/AMPManager.h"
#include "utils/InputManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>
#include "materials/Material.h"

#include "utils/AMPManager.h"
#include "utils/ReadTestMesh.h"
#include "utils/WriteSolutionToFile.h"

#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshVariable.h"
#include "ampmesh/MeshAdapter.h"

#include "operators/ColumnOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/mechanics/MechanicsLinearFEOperator.h"

#include "mesh_communication.h"

#include "MPCtestUtils.h"

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);
  std::string mesh_file1 = input_db->getString("mesh_file1");

  AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );

  const unsigned int mesh_dim = 3;
  boost::shared_ptr< ::Mesh > meshMaster(new ::Mesh(mesh_dim));
  boost::shared_ptr< ::Mesh > meshSlave(new ::Mesh(mesh_dim));

  bool binaryMeshes = input_db->getBool("BinaryMeshes");

  if(binaryMeshes) {
    AMP::readBinaryTestMesh(mesh_file1, meshMaster);
    AMP::readBinaryTestMesh(mesh_file1, meshSlave);
  } else {
    AMP::readTestMesh(mesh_file1, meshMaster);
    AMP::readTestMesh(mesh_file1, meshSlave);
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
										   "BVPOperator",
										   input_db,
										   masterElementPhysicsModel));

  AMP::LinearAlgebra::Matrix::shared_ptr masterMatrix = masterOperator->getMatrix();

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel;
  boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> slaveOperator = boost::dynamic_pointer_cast<
  AMP::Operator::MechanicsLinearFEOperator>(AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter,
											   "VolumeOperator",
											   input_db,
											   slaveElementPhysicsModel));

  AMP::LinearAlgebra::Matrix::shared_ptr slaveMatrix = slaveOperator->getMatrix();

  boost::shared_ptr<AMP::Operator::OperatorParameters> dummyParams;
  boost::shared_ptr<AMP::Operator::ColumnOperator> columnOperator(new AMP::Operator::ColumnOperator(dummyParams));
  columnOperator->append(masterOperator);
  columnOperator->append(slaveOperator);

  AMP::LinearAlgebra::Variable::shared_ptr masterVar = masterOperator->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr slaveVar = slaveOperator->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr columnVar = columnOperator->getOutputVariable();

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOperator1 =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
									  AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter, "LoadOperator", input_db,  dummyModel));
  loadOperator1->setVariable(slaveVar);

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr columnSolVec = manager->createVector(columnVar);
  AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec = columnSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr columnResVec = columnSolVec->cloneVector();

  AMP::LinearAlgebra::Vector::shared_ptr masterSolVec = columnSolVec->subsetVectorForVariable(masterVar);
  AMP::LinearAlgebra::Vector::shared_ptr slaveSolVec = columnSolVec->subsetVectorForVariable(slaveVar);

  columnSolVec->zero();
  columnRhsVec->zero();
  columnResVec->zero();

  loadOperator1->apply(nullVec, nullVec, columnRhsVec, 1.0, 0.0);

  int maxIters = input_db->getInteger("maxIters");

  {
    std::cout<<"MPC-PCG: "<<std::endl;
    AMP::LinearAlgebra::Vector::shared_ptr MatOutVec = columnSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr solOldVec = columnSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr resOldVec = columnSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr pOldVec = columnSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr pNewVec = columnSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr zVec = columnSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr diagMaster = masterMatrix->extractDiagonal();
    AMP::LinearAlgebra::Vector::shared_ptr diagSlave = slaveMatrix->extractDiagonal();
    AMP::LinearAlgebra::Vector::shared_ptr diagMasterInverse = diagMaster->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr diagSlaveInverse = diagSlave->cloneVector();

    diagMasterInverse->reciprocal(diagMaster);
    diagSlaveInverse->reciprocal(diagSlave);

    columnOperator->apply(nullVec, columnSolVec, MatOutVec, 1.0, 0.0);
    addSlaveToMaster(slaveVar, masterVar, MatOutVec, slaveMeshAdapter, masterMeshAdapter, slaveNodes, masterNodes);
    setSlaveToZero(slaveVar, MatOutVec, slaveMeshAdapter, slaveNodes);

    columnResVec->subtract(columnRhsVec, MatOutVec);

    AMP::LinearAlgebra::Vector::shared_ptr zMasterVec = zVec->subsetVectorForVariable(masterVar);
    AMP::LinearAlgebra::Vector::shared_ptr zSlaveVec = zVec->subsetVectorForVariable(slaveVar);
    AMP::LinearAlgebra::Vector::shared_ptr resMasterVec = columnResVec->subsetVectorForVariable(masterVar);
    AMP::LinearAlgebra::Vector::shared_ptr resSlaveVec = columnResVec->subsetVectorForVariable(slaveVar);
    zMasterVec->multiply(diagMasterInverse, resMasterVec);
    zSlaveVec->multiply(diagSlaveInverse, resSlaveVec);

    pOldVec->copyVector(zVec);

    for(int iter = 0; iter < maxIters; iter++) {
      double resNorm = columnResVec->L2Norm();
      double solMaxNorm = columnSolVec->maxNorm();

      setSlaveToZero(slaveVar, columnSolVec, slaveMeshAdapter, slaveNodes);
      double solNorm2 = columnSolVec->L2Norm();
      std::cout<<"Iter = "<<iter<<" ResNorm2 = "<<std::setprecision(15)<<resNorm<<
        " SolMaxNorm = "<<std::setprecision(15)<<solMaxNorm<<
        " SolNorm2 = "<<std::setprecision(15)<<solNorm2<<std::endl;
      copyMasterToSlave(slaveVar, masterVar, columnSolVec, slaveMeshAdapter, masterMeshAdapter, slaveNodes, masterNodes);

      copyMasterToSlave(slaveVar, masterVar, pOldVec, slaveMeshAdapter, masterMeshAdapter, slaveNodes, masterNodes);

      resOldVec->copyVector(columnResVec);
      solOldVec->copyVector(columnSolVec);

      double resOldDotZ = resOldVec->dot(zVec);

      columnOperator->apply(nullVec, pOldVec, MatOutVec, 1.0, 0.0);
      addSlaveToMaster(slaveVar, masterVar, MatOutVec, slaveMeshAdapter, masterMeshAdapter, slaveNodes, masterNodes);
      setSlaveToZero(slaveVar, MatOutVec, slaveMeshAdapter, slaveNodes);

      double alphaDenom = MatOutVec->dot(pOldVec);

      double alpha = resOldDotZ/alphaDenom;

      columnSolVec->axpy(alpha, pOldVec, solOldVec);

      columnResVec->axpy(-alpha, MatOutVec, resOldVec);

      zMasterVec = zVec->subsetVectorForVariable(masterVar);
      zSlaveVec = zVec->subsetVectorForVariable(slaveVar);
      resMasterVec = columnResVec->subsetVectorForVariable(masterVar);
      resSlaveVec = columnResVec->subsetVectorForVariable(slaveVar);
      zMasterVec->multiply(diagMasterInverse, resMasterVec);
      zSlaveVec->multiply(diagSlaveInverse, resSlaveVec);

      double resNewDotZ = columnResVec->dot(zVec);

      double beta = resNewDotZ/resOldDotZ;

      std::cout<<"Iter = "<<iter<<" alpha = "<<std::setprecision(15)<<alpha<<
        " beta = "<<std::setprecision(15)<<beta<<std::endl<<std::endl;

      pNewVec->axpy(beta, pOldVec, zVec);

      pOldVec->copyVector(pNewVec);
    }
  }

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
    exeNames.push_back("testPCG_MPC-1");
    exeNames.push_back("testPCG_MPC-2");
    exeNames.push_back("testPCG_MPC-3");
    exeNames.push_back("testPCG_MPC-4");
    exeNames.push_back("testPCG_MPC-5");
  } else {
    for(int i = 1; i < argc; i++) {
      char inpName[100];
      sprintf(inpName, "testPCG_MPC-%s", argv[i]);
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

