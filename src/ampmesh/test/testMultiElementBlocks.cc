#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <string>
#include "utils/AMPManager.h"
#include "ampmesh/MeshManager.h"
#include "boost/shared_ptr.hpp"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshAdapter.h"

#include "ampmesh/SiloIO.h"

#include "elem.h"
#include <exception>

void multiBlockTest(AMP::UnitTest *ut , std::string exeName)
{
  // Initialization
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logAllNodes(log_file);

  std::cout << "testing with input file " << input_file << std::endl;
  std::cout.flush();

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ( "multiblock" );

  AMP::Mesh::MeshManager::Adapter::ElementIterator elem = meshAdapter->beginElement();
  AMP::Mesh::MeshManager::Adapter::ElementIterator end_elem = meshAdapter->endElement();
  
  for( ; elem != end_elem ; ++elem)
  {
      int elemBlock;
      ::Elem * elemPtr = &(elem->getElem()) ; 

      elemBlock = elemPtr->subdomain_id() ;

      std::cout << "Element Block ID is : " << elemBlock <<std::endl;
  }

}

  
// Input and output file names
int main(int argc, char *argv[])
{

    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);

    AMP::UnitTest ut;

    const int NUMFILES=1;
    std::string files[NUMFILES] = {
        "testMultiElementBlock-1"
    };

    try {
        for (int i=0; i<NUMFILES; i++) {
            multiBlockTest(&ut, files[i]);
        }
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


