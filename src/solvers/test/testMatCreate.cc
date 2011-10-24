
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

void myTest(AMP::UnitTest *ut, std::string input_file) {
  std::string log_file = "output_testMatCreate";

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);

  AMP::Mesh::MeshManagerParameters::shared_ptr meshmgrParams (new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager (new AMP::Mesh::MeshManager ( meshmgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ("cube");

  AMP::LinearAlgebra::Variable::shared_ptr inVar (new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3>("inputVar", meshAdapter));
  AMP::LinearAlgebra::Variable::shared_ptr outVar (new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>("outputVar", meshAdapter));

  AMP::LinearAlgebra::Matrix::shared_ptr mat1 = meshAdapter->createMatrix(inVar, outVar);

  if(mat1.get()!=NULL){
    ut->passes("Able to create a non-square matrices");
  }else{
    ut->failure("Unable to create a non-square matrices");
  }

}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  try {
    myTest(&ut, "input_testLinearFlow-1");
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



