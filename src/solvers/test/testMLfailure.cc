
#include <iostream>

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "mesh_communication.h"

#include "mpi.h"

#include "boost/shared_ptr.hpp"

#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>


#include "operators/OperatorBuilder.h"
#include "operators/LinearBVPOperator.h"
#include "operators/boundary/DirichletVectorCorrection.h"

#include "solvers/PetscKrylovSolverParameters.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/TrilinosMLSolver.h"

#include "utils/ReadTestMesh.h"

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  boost::shared_ptr<AMP::Database> ml_db = input_db->getDatabase("ML_Solver"); 
  boost::shared_ptr<AMP::Database> cg_db = input_db->getDatabase("CG_Solver"); 

  AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );

  std::string meshFile = input_db->getString("mesh");

  const unsigned int mesh_dim = 3;
  boost::shared_ptr< ::Mesh > mesh(new ::Mesh(mesh_dim));

  if(globalComm.getRank() == 0) {
    AMP::readTestMesh(meshFile, mesh);
  }

  MeshCommunication().broadcast(*(mesh.get()));
  mesh->prepare_for_use(false);

  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter ( new AMP::Mesh::MeshManager::Adapter (mesh) );

  manager->addMesh(meshAdapter, "mesh");

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator = boost::dynamic_pointer_cast<
    AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "BVPOperator", input_db, elementPhysicsModel));

  AMP::LinearAlgebra::Variable::shared_ptr dispVar = bvpOperator->getOutputVariable();

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOperator = boost::dynamic_pointer_cast<
    AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "LoadOperator", input_db, dummyModel));
  loadOperator->setVariable(dispVar);

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr solVec = manager->createVector(dispVar);
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = manager->createVector(dispVar);
  AMP::LinearAlgebra::Vector::shared_ptr resVec = manager->createVector(dispVar);

  rhsVec->zero();
  loadOperator->apply(nullVec, nullVec, rhsVec, 1.0, 0.0);

  solVec->zero();
  resVec->zero();

  if(globalComm.getRank() == 0) {
    std::cout<<"Solving using ML"<<std::endl;
  }

  boost::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> mlParams(new AMP::Solver::TrilinosMLSolverParameters(ml_db));
  mlParams->d_pOperator = bvpOperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> mlPC(new AMP::Solver::TrilinosMLSolver(mlParams));

  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> cgParams(new AMP::Solver::PetscKrylovSolverParameters(cg_db));
  cgParams->d_pOperator = bvpOperator;
  cgParams->d_comm = globalComm;
  cgParams->d_pPreconditioner = mlPC;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> cgSolver(new AMP::Solver::PetscKrylovSolver(cgParams));

  cgSolver->solve(rhsVec, solVec);

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testMLfailure";

  try {
    myTest(&ut, exeName);
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


