
#include <iostream>
#include <string>

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/ReadTestMesh.h"

#include "ampmesh/MeshVariable.h"
#include "ampmesh/MeshAdapter.h"
#include "ampmesh/MeshManager.h"

#include "mesh_communication.h"

void myTest(AMP::UnitTest *ut, std::string mesh_file) {
  std::string log_file = "output_testMatMultiply";

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> mesh_file_db(new AMP::InputDatabase("mesh_file_db"));
  AMP::InputManager::getManager()->parseInputFile(mesh_file, mesh_file_db);

  const unsigned int mesh_dim = 3;
  boost::shared_ptr< ::Mesh > myMesh(new ::Mesh(mesh_dim));

  AMP::readTestMesh(mesh_file, myMesh);

  MeshCommunication().broadcast(*(myMesh.get()));

  myMesh->prepare_for_use(false);

  AMP::Mesh::MeshManager::Adapter::shared_ptr myMeshAdapter ( new AMP::Mesh::MeshManager::Adapter (myMesh) );

  AMP::LinearAlgebra::Variable::shared_ptr myVar (new 
      AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3>("myVar", myMeshAdapter));

  AMP::LinearAlgebra::Vector::shared_ptr vec1 = myMeshAdapter->createVector(myVar);
  vec1->setToScalar(1.0);

  AMP::LinearAlgebra::Matrix::shared_ptr mat1 = myMeshAdapter->createMatrix(myVar, myVar);
  mat1->zero();
  mat1->setDiagonal(vec1);

  AMP::LinearAlgebra::Matrix::shared_ptr mat2 = mat1->cloneMatrix();
  mat2->zero();
  mat2->setDiagonal(vec1);

  AMP::LinearAlgebra::Matrix::shared_ptr mat3 = AMP::LinearAlgebra::Matrix::matMultiply(mat1, mat2);

  std::vector<unsigned int> cols1;
  std::vector<double> vals1;
  mat1->getRowByGlobalID(0, cols1, vals1);

  std::vector<unsigned int> cols2;
  std::vector<double> vals2;
  mat2->getRowByGlobalID(0, cols2, vals2);

  std::vector<unsigned int> cols3;
  std::vector<double> vals3;
  mat3->getRowByGlobalID(0, cols3, vals3);

  ut->passes("testMatMultiply");

}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  try {
    myTest(&ut, "mpcMesh-1");
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



