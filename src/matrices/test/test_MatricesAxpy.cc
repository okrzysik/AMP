#include <iostream>

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "ampmesh/Mesh.h"

#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"
#include "matrices/MatrixBuilder.h"


void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;
  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  // build mesh database 
  AMP::shared_ptr<AMP::Database> database(new AMP::MemoryDatabase("Mesh"));
  database->putString("MeshName", "cube");
  database->putString("MeshType", "AMP");
  database->putInteger("dim", 3);
  database->putString("Generator", "cube");
  database->putDouble("x_offset", 0.0);
  database->putDouble("y_offset", 0.0);
  database->putDouble("z_offset", 0.0);
  std::vector<int> size(3, 8);
  database->putIntegerArray("Size", size);
  std::vector<double> range(6, 0.0);
  range[1] = range[3] = range[5] = 1.0;
  database->putDoubleArray("Range", range);
  // create mesh
  AMP::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
  params->setComm(globalComm);
  AMP::shared_ptr<AMP::Mesh::Mesh> mesh = AMP::Mesh::Mesh::buildMesh(params); 
  // create two different dof managers
  AMP::Discretization::DOFManager::shared_ptr firstDofManager = AMP::Discretization::simpleDOFManager::create(mesh, AMP::Mesh::Vertex, 1, 3, true);
  AMP::Discretization::DOFManager::shared_ptr secondDofManager = AMP::Discretization::simpleDOFManager::create(mesh, AMP::Mesh::Volume, 0, 8, true);
  // create the two corresponding vectors
  AMP::LinearAlgebra::Variable::shared_ptr var(new  AMP::LinearAlgebra::Variable("var"));
  AMP::LinearAlgebra::Vector::shared_ptr firstVec = AMP::LinearAlgebra::createVector(firstDofManager, var, true); // n
  AMP::LinearAlgebra::Vector::shared_ptr secondVec = AMP::LinearAlgebra::createVector(secondDofManager, var, true); // m
  firstVec->zero();
  secondVec->zero();
  // create four matrices
  AMP::LinearAlgebra::Matrix::shared_ptr firstMat = AMP::LinearAlgebra::createMatrix(firstVec, firstVec); // nxn
  AMP::LinearAlgebra::Matrix::shared_ptr secondMat = AMP::LinearAlgebra::createMatrix(secondVec, secondVec); // mxm
  firstMat->setScalar(1.0);
  secondMat->setScalar(2.0);
  // axpy
  try {
    firstMat->axpy(-1.0, secondMat);
    std::cerr<<"this should not be allowed..."<<std::endl;
    ut->failure(exeName);
  } catch (std::exception & err) {
    ut->expected_failure(exeName);
  } catch (...) {
    std::cerr<<"unknown error"<<std::endl;
    ut->failure(exeName);
  }
 
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "test_MatricesAxpy";

  try {
    myTest(&ut, exeName);
//    ut.passes(exeName);
  } catch (std::exception &err) {
    std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
    ut.failure("ERROR: While testing");
  } catch( ... ) {
    std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
    ut.failure("ERROR: While testing");
  }

  ut.report();

  int numFailed = ut.NumFailGlobal();
  AMP::AMPManager::shutdown();
  return numFailed;
}
