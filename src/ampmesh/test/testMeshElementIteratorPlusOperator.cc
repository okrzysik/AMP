
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"

void myTest(AMP::UnitTest *ut) {
  std::string input_file = "input_Mesh";
  std::string log_file = "output_testMeshElementIteratorPlusOperator"; 

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  // Load the input file
  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  // Load the mesh
  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);

  AMP::Mesh::MeshIterator el = meshAdapter->getIterator(AMP::Mesh::Volume, 0);
  AMP::Mesh::MeshIterator end_el = el.end();
  AMP_ASSERT(el != end_el);
  for(; el != end_el; ++el) {
    std::vector<AMP::Mesh::MeshElement> vertices = el->getElements(AMP::Mesh::Vertex);
  }//end el
  ut->passes("testMeshElementIteratorPlusOperator-1a");
  std::cout<<"Passed 1A"<<std::endl;

  size_t localNumElems = meshAdapter->numLocalElements(AMP::Mesh::Volume);

  el = meshAdapter->getIterator(AMP::Mesh::Volume, 0);
  for(size_t i = 0; i < localNumElems; ++i) {
    std::vector<AMP::Mesh::MeshElement> vertices = (el + i)->getElements(AMP::Mesh::Vertex);
  }//end i
  ut->passes("testMeshElementIteratorPlusOperator-1b");
  std::cout<<"Passed 1B"<<std::endl;

  el = meshAdapter->getIterator(AMP::Mesh::Volume, 0);
  end_el = el.end();
  for(; el != end_el; ++el) {
    std::vector<AMP::Mesh::MeshElement> vertices = (&(*el))->getElements(AMP::Mesh::Vertex);
  }//end el
  ut->passes("testMeshElementIteratorPlusOperator-2a");
  std::cout<<"Passed 2A"<<std::endl;

  el = meshAdapter->getIterator(AMP::Mesh::Volume, 0);
  for(size_t i = 0; i < localNumElems; ++i) {
    std::vector<AMP::Mesh::MeshElement> vertices = (&(*(el + i)))->getElements(AMP::Mesh::Vertex);
  }//end i
  ut->passes("testMeshElementIteratorPlusOperator-2b");
  std::cout<<"Passed 2B"<<std::endl;

  el = meshAdapter->getIterator(AMP::Mesh::Volume, 0);
  end_el = el.end();
  for(; el != end_el; ++el) {
    AMP::Mesh::MeshElement* elem = (&(*el));
    std::vector<AMP::Mesh::MeshElement> vertices = elem->getElements(AMP::Mesh::Vertex);
  }//end el
  ut->passes("testMeshElementIteratorPlusOperator-3a");
  std::cout<<"Passed 3A"<<std::endl;

  el = meshAdapter->getIterator(AMP::Mesh::Volume, 0);
  for(size_t i = 0; i < localNumElems; ++i) {
    AMP::Mesh::MeshElement* elem = (&(*(el + i)));
    std::vector<AMP::Mesh::MeshElement> vertices = elem->getElements(AMP::Mesh::Vertex);
  }//end i
  ut->passes("testMeshElementIteratorPlusOperator-3b");
  std::cout<<"Passed 3B"<<std::endl;

}


int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  try {
    myTest(&ut);
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


