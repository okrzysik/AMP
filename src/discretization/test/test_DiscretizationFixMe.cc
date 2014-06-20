#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "ampmesh/Mesh.h"

#include "discretization/simpleDOF_Manager.h"

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;
  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  boost::shared_ptr<AMP::Database> database(new AMP::MemoryDatabase("Mesh"));
  // build mesh database
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
  boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
  params->setComm(globalComm);
  boost::shared_ptr<AMP::Mesh::Mesh> mesh = AMP::Mesh::Mesh::buildMesh(params);
  // create dof manager
  AMP::Discretization::DOFManager::shared_ptr dofManager = AMP::Discretization::simpleDOFManager::create(mesh, AMP::Mesh::Vertex, 1, 1, true);
  AMP::Mesh::MeshIterator meshIterator = mesh->getIterator(AMP::Mesh::Volume, 0);
  size_t const numLocalElements = mesh->numLocalElements(AMP::Mesh::Volume);
  std::vector<size_t> dofIndices;
  for (size_t elem = 0; elem < numLocalElements; ++elem, ++meshIterator) {
    // we have to do the following to get the dof indices
    std::vector<AMP::Mesh::MeshElement> verticesInElement = meshIterator->getElements(AMP::Mesh::Vertex);
    size_t const numVertices = verticesInElement.size();
    std::vector<AMP::Mesh::MeshElementID> verticesIDs(numVertices);
    for (size_t vertex = 0; vertex < numVertices; ++vertex) {
      verticesIDs[vertex] = verticesInElement[vertex].globalID();
    } // end for point
    dofManager->getDOFs(verticesIDs, dofIndices);
    size_t const numDOFsInElement = dofIndices.size();
    // this is silly
    // we should be able to do
    dofManager->getDOFs(meshIterator->globalID(), dofIndices);
    AMP_ASSERT(dofIndices.size() == numDOFsInElement);
  } // end for elem

  ut->passes(exeName);
}



int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testDamien";

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
