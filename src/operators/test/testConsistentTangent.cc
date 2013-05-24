
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
#include <cstdlib>

#include "ampmesh/Mesh.h"
#include "ampmesh/libmesh/libMesh.h"

#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "libmesh.h"
#include "mesh_communication.h"

#include "operators/OperatorBuilder.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"

#include "utils/ReadTestMesh.h"

void myTest(AMP::UnitTest *ut, std::string exeName, int callLinReset) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;
  std::string msgName = exeName;
  if(callLinReset) {
    log_file = log_file + "-1";
    msgName = msgName + "-1";
  } else {
    log_file = log_file + "-0";
    msgName = msgName + "-0";
  }

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  const unsigned int mesh_dim = 3;
  boost::shared_ptr< ::Mesh > mesh(new ::Mesh(mesh_dim));

  std::string mesh_file = input_db->getString("mesh_file");

  if( globalComm.getRank() == 0 ) {
    AMP::readTestMesh(mesh_file, mesh);
  }//end if root processor

  MeshCommunication().broadcast(*(mesh.get()));

  mesh->prepare_for_use(false);

  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::shared_ptr (
      new AMP::Mesh::libMesh (mesh, "TestMesh") );

  boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> nonlinOperator =
    boost::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
        AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "NonlinearMechanicsOperator", input_db));
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel = nonlinOperator->getMaterialModel();

  boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> linOperator =
    boost::dynamic_pointer_cast<AMP::Operator::MechanicsLinearFEOperator>(
        AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "LinearMechanicsOperator", input_db, elementPhysicsModel));


  AMP::Discretization::DOFManager::shared_ptr dofMap = AMP::Discretization::simpleDOFManager::create(
      meshAdapter, AMP::Mesh::Vertex, 1, 3, true); 

  AMP::LinearAlgebra::Variable::shared_ptr var = nonlinOperator->getOutputVariable();

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector(dofMap, var, true);
  AMP::LinearAlgebra::Vector::shared_ptr resVecNonlin = solVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr resVecLin = solVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr resDiffVec = solVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr tmpNonlinVec = solVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr tmpLinVec = solVec->cloneVector();

  solVec->setToScalar(0.0);
  double solNorm = solVec->L2Norm();
  AMP::pout<<"sol-Norm-2 = "<<solNorm<<std::endl;

  nonlinOperator->apply(nullVec, solVec, resVecNonlin, 1.0, 0.0);
  linOperator->reset(nonlinOperator->getJacobianParameters(solVec));
  linOperator->apply(nullVec, solVec, resVecLin, 1.0, 0.0);
  resDiffVec->subtract(resVecNonlin, resVecLin);

  double epsilon = 1.0e-13*(((linOperator->getMatrix())->extractDiagonal())->L1Norm());
  AMP::pout<<"epsilon = "<<epsilon<<std::endl;

  double nonLinNorm = resVecNonlin->L1Norm();
  double linNorm = resVecLin->L1Norm();
  double diffNorm = resDiffVec->L1Norm();

  AMP::pout<<"nonLin-Norm-1 = "<<nonLinNorm<<std::endl;
  AMP::pout<<"lin-Norm-1 = "<<linNorm<<std::endl;
  AMP::pout<<"diff-Norm-1 = "<<diffNorm<<std::endl;

  if( nonLinNorm > epsilon ) {
    ut->failure(msgName);
  }

  if( linNorm > epsilon ) {
    ut->failure(msgName);
  }

  if( diffNorm > epsilon ) {
    ut->failure(msgName);
  }

  solVec->setRandomValues();
  solVec->scale(100.0);
  solVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
  solNorm = solVec->L2Norm();
  AMP::pout<<"sol-Norm-2 = "<<solNorm<<std::endl;

  nonlinOperator->apply(nullVec, solVec, resVecNonlin, 1.0, 0.0);
  if(callLinReset) {
    linOperator->reset(nonlinOperator->getJacobianParameters(solVec));
  }
  linOperator->apply(nullVec, solVec, resVecLin, 1.0, 0.0);
  resDiffVec->subtract(resVecNonlin, resVecLin);

  nonLinNorm = resVecNonlin->L2Norm();
  linNorm = resVecLin->L2Norm();
  diffNorm = resDiffVec->L1Norm();

  AMP::pout<<"nonLin-Norm-2 = "<<nonLinNorm<<std::endl;
  AMP::pout<<"lin-Norm-2 = "<<linNorm<<std::endl;
  AMP::pout<<"diff-Norm-1 = "<<diffNorm<<std::endl;

  if( diffNorm > epsilon ) {
    ut->failure(msgName);
  }

  tmpNonlinVec->copyVector(resVecNonlin);
  tmpNonlinVec->scale(10.0);

  tmpLinVec->copyVector(resVecLin);
  tmpLinVec->scale(10.0);

  solVec->scale(10.0);
  solVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
  solNorm = solVec->L2Norm();
  AMP::pout<<"sol-Norm-2 = "<<solNorm<<std::endl;

  nonlinOperator->apply(nullVec, solVec, resVecNonlin, 1.0, 0.0);
  if(callLinReset) {
    linOperator->reset(nonlinOperator->getJacobianParameters(solVec));
  }
  linOperator->apply(nullVec, solVec, resVecLin, 1.0, 0.0);

  nonLinNorm = resVecNonlin->L2Norm();
  linNorm = resVecLin->L2Norm();
  AMP::pout<<"nonLin-Norm-2 = "<<nonLinNorm<<std::endl;
  AMP::pout<<"lin-Norm-2 = "<<linNorm<<std::endl;

  resDiffVec->subtract(resVecNonlin, tmpNonlinVec);
  diffNorm = resDiffVec->L1Norm();

  if( diffNorm > (10.0*epsilon) ) {
    ut->failure(msgName);
  }

  resDiffVec->subtract(resVecLin, tmpLinVec);
  diffNorm = resDiffVec->L1Norm();

  if( diffNorm > (10.0*epsilon) ) {
    ut->failure(msgName);
  }

  ut->passes(msgName);
}

int main(int argc, char *argv[]) {

  AMP::AMPManager::startup(argc, argv);
  boost::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit(new AMP::Mesh::initializeLibMesh(AMP_COMM_WORLD));

  AMP::UnitTest ut;

  std::vector<std::string> exeNames;

  exeNames.push_back("testConsistentTangent-1-mesh1-normal");
  exeNames.push_back("testConsistentTangent-1-mesh1-reduced");
  exeNames.push_back("testConsistentTangent-1-cookMesh1-normal");
  exeNames.push_back("testConsistentTangent-1-cookMesh1-reduced");

  exeNames.push_back("testConsistentTangent-2-mesh1-normal");
  exeNames.push_back("testConsistentTangent-2-mesh1-reduced");
  exeNames.push_back("testConsistentTangent-2-cookMesh1-normal");
  exeNames.push_back("testConsistentTangent-2-cookMesh1-reduced");

  exeNames.push_back("testConsistentTangent-3-mesh1-normal");
  exeNames.push_back("testConsistentTangent-3-mesh1-reduced");
  exeNames.push_back("testConsistentTangent-3-cookMesh1-normal");
  exeNames.push_back("testConsistentTangent-3-cookMesh1-reduced");

  exeNames.push_back("testConsistentTangent-4-mesh1-normal");
  exeNames.push_back("testConsistentTangent-4-mesh1-reduced");
  exeNames.push_back("testConsistentTangent-4-cookMesh1-normal");
  exeNames.push_back("testConsistentTangent-4-cookMesh1-reduced");

  for(int j = 0; j < 2; j++ ) {
    for(size_t i = 0; i < exeNames.size(); i++) {
      try {
        myTest(&ut, exeNames[i], j);
      } catch (std::exception &err) {
        AMP::pout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
        ut.failure("ERROR: While testing");
      } catch( ... ) {
        AMP::pout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
        ut.failure("ERROR: While testing");
      }
    }//end for i
  }

  ut.report();
  int num_failed = ut.NumFailGlobal();

  libmeshInit.reset();
  AMP::AMPManager::shutdown();
  return num_failed;
}


