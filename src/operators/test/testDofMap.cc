
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
#include <vector>

#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshVariable.h"


void myTest(AMP::UnitTest *ut, std::string exeName)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST( input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!" );
  std::string mesh_file = input_db->getString("Mesh");

  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter =
    AMP::Mesh::MeshManager::Adapter::shared_ptr(new AMP::Mesh::MeshManager::Adapter());
  meshAdapter->readExodusIIFile ( mesh_file.c_str() );

  AMP::LinearAlgebra::Variable::shared_ptr var(new AMP::LinearAlgebra::VectorVariable<
      AMP::Mesh::NodalVariable, 1>("temp", meshAdapter)); 

  AMP::Mesh::DOFMap::shared_ptr dofMap = meshAdapter->getDOFMap(var);

  AMP::Mesh::MeshManager::Adapter::ElementIterator  el = meshAdapter->beginElement();
  AMP::Mesh::MeshManager::Adapter::ElementIterator  end_el = meshAdapter->endElement();

  for(int elId = 0; el != end_el; ++el, ++elId) {
    std::vector<unsigned int> dofIds; 
    dofMap->getDOFs(*el, dofIds, 0);
    AMP_ASSERT(dofIds.size() == 8);
    std::cout<<"e = "<<elId<<" : "<<std::endl;
    for(int i = 0; i < 8; ++i) {
      AMP::Mesh::LibMeshPoint pt = el->getPoint(i);
      std::cout<<"v = "<<i<<" : ";
      std::cout<<"x = "<<(pt.x())<<" : ";
      std::cout<<"y = "<<(pt.y())<<" : ";
      std::cout<<"z = "<<(pt.z())<<" : ";
      std::cout<<"d0 = "<<(dofIds[i])<<std::endl;
    }//end i
    std::cout<<std::endl;
  }//end el

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::vector<std::string> exeNames;
  exeNames.push_back("testDofMap");

  for(unsigned int i = 0; i < exeNames.size(); i++) {
    try {
      myTest(&ut, exeNames[i]);
    } catch (std::exception &err) {
      std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
      ut.failure("ERROR: While testing");
    } catch( ... ) {
      std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
      ut.failure("ERROR: While testing");
    }
  }

  ut.report();

  int num_failed = ut.NumFailGlobal();
  AMP::AMPManager::shutdown();
  return num_failed;

}   


