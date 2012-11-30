
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
  std::string input_file = "input_testGetElements";
  std::string log_file = "output_testGetElements"; 

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  const int rank = globalComm.getRank();

  // Parse the input file
  globalComm.barrier();
  double time1b = MPI_Wtime();
  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);
  globalComm.barrier();
  double time1e = MPI_Wtime();
  if(!rank) {
    std::cout<<"Time to parse input file = "<<(time1e - time1b)<<std::endl;
  }

  // Load the mesh
  globalComm.barrier();
  double time2b = MPI_Wtime();
  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);
  globalComm.barrier();
  double time2e = MPI_Wtime();
  if(!rank) {
    std::cout<<"Time to load the mesh = "<<(time2e - time2b)<<std::endl;
  }

  size_t localNumElems = meshAdapter->numLocalElements(AMP::Mesh::Volume);

  globalComm.barrier();
  double time3b = MPI_Wtime();
  AMP::Mesh::MeshIterator el = meshAdapter->getIterator(AMP::Mesh::Volume, 0);
  for(size_t i = 0; i < localNumElems; ++i, ++el) {
    std::vector<AMP::Mesh::MeshElement> vertices = el->getElements(AMP::Mesh::Vertex);
  }//end 
  globalComm.barrier();
  double time3e = MPI_Wtime();
  if(!rank) {
    std::cout<<"Time to iterate (Type 1) through the mesh elements and access the vertices = "<<(time3e - time3b)<<std::endl;
  }

  globalComm.barrier();
  double time4b = MPI_Wtime();
  el = meshAdapter->getIterator(AMP::Mesh::Volume, 0);
  for(size_t i = 0; i < localNumElems; ++i) {
    std::vector<AMP::Mesh::MeshElement> vertices = (el + i)->getElements(AMP::Mesh::Vertex);
  }//end 
  globalComm.barrier();
  double time4e = MPI_Wtime();
  if(!rank) {
    std::cout<<"Time to iterate (Type 2) through the mesh elements and access the vertices = "<<(time4e - time4b)<<std::endl;
  }


  ut->passes("testGetElements");
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



