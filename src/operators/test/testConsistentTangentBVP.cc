#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>
#include <cstdlib>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "vectors/PetscVector.h"

#include "MeshManager.h"
#include "MeshVariable.h"

#include "libmesh.h"
#include "mesh_communication.h"

#include "materials/Material.h"
#include "OperatorBuilder.h"
#include "LinearBVPOperator.h"
#include "NonlinearBVPOperator.h"

#include "ReadTestMesh.h"


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
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );

  std::string mesh_file = input_db->getString("mesh_file");

  const unsigned int mesh_dim = 3;
  boost::shared_ptr< ::Mesh > mesh(new ::Mesh(mesh_dim));

  if(ut->rank() == 0) {
    AMP::readTestMesh(mesh_file, mesh);
  }//end if root processor

  MeshCommunication().broadcast(*(mesh.get()));

  mesh->prepare_for_use(false);

  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter ( new AMP::Mesh::MeshManager::Adapter (mesh) );

  manager->addMesh(meshAdapter, "mesh");

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinOperator =
    boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
														    "NonlinearMechanicsOperator",
														    input_db,														    
														    elementPhysicsModel));

  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linOperator =
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
														 "LinearMechanicsOperator",
														 input_db,
														 elementPhysicsModel));

  AMP::LinearAlgebra::Variable::shared_ptr dispVar = nonlinOperator->getOutputVariable();

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr solVec = meshAdapter->createVector(dispVar);
  AMP::LinearAlgebra::Vector::shared_ptr resVecNonlin = meshAdapter->createVector(dispVar);
  AMP::LinearAlgebra::Vector::shared_ptr resVecLin = meshAdapter->createVector(dispVar);
  AMP::LinearAlgebra::Vector::shared_ptr resDiffVec = meshAdapter->createVector(dispVar);
  AMP::LinearAlgebra::Vector::shared_ptr tmpNonlinVec = meshAdapter->createVector(dispVar);
  AMP::LinearAlgebra::Vector::shared_ptr tmpLinVec = meshAdapter->createVector(dispVar);

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
  nonlinOperator->modifyInitialSolutionVector(solVec);
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
  nonlinOperator->modifyInitialSolutionVector(solVec);
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
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;

    exeNames.push_back("testConsistentTangentBVP-1-mesh1-normal");
    exeNames.push_back("testConsistentTangentBVP-2-mesh1-normal");
    exeNames.push_back("testConsistentTangentBVP-3-mesh1-normal");
    exeNames.push_back("testConsistentTangentBVP-4-mesh1-normal");

    exeNames.push_back("testConsistentTangentBVP-4-mesh1-reduced");
  
    for(int j = 0; j < 2; j++ ) {
        for(size_t i = 0; i < exeNames.size(); i++) {
            try {
                myTest(&ut, exeNames[i], j);
            } catch (std::exception &err) {
                AMP::pout << "ERROR: While testing "<<argv[0]<< err.what() << std::endl;
                ut.failure("ERROR: While testing");
            } catch( ... ) {
                AMP::pout << "ERROR: While testing "<<argv[0]<< "An unknown exception was thrown." << std::endl;
                ut.failure("ERROR: While testing");
            }
        }//end for i
    }


    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}


