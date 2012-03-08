
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshVariable.h"

void myTest(AMP::UnitTest *ut)
{
  std::string exeName("testGhostWrite");
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  std::string mesh_file = input_db->getString("Mesh");

  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = AMP::Mesh::MeshManager::Adapter::shared_ptr ( new AMP::Mesh::MeshManager::Adapter () );
  meshAdapter->readExodusIIFile ( mesh_file.c_str() );

  AMP::LinearAlgebra::Variable::shared_ptr var(new AMP::Mesh::NodalScalarVariable("dummy", meshAdapter)); 
  AMP::LinearAlgebra::Matrix::shared_ptr mat = meshAdapter->getMatrix(var, var); 

  AMP::Mesh::DOFMap::shared_ptr dof_map = meshAdapter->getDOFMap ( var );

  AMP::Mesh::MeshManager::Adapter::NodeIterator nd = meshAdapter->beginNode(  );
  AMP::Mesh::MeshManager::Adapter::NodeIterator end_nd = meshAdapter->endNode(  );

  for( ; nd != end_nd; ++nd) {
    AMP::Mesh::MeshManager::Adapter::NodeElementIterator el =  meshAdapter->beginElementForNode ( *nd );
    AMP::Mesh::MeshManager::Adapter::NodeElementIterator end_el = meshAdapter->endElementForNode ( *nd );

    std::vector<unsigned int> ndGlobalIds;
    std::vector<unsigned int> emptyVec;
    dof_map->getDOFs(*nd, ndGlobalIds, emptyVec);

    for( ; el != end_el; ++el) {
      std::vector<unsigned int> dofIndices;
      dof_map->getDOFs(*el, dofIndices);
      for(unsigned int j = 0; j < ndGlobalIds.size(); j++) {
        for(unsigned int i = 0; i < dofIndices.size(); i++) {
          mat->setValueByGlobalID ( ndGlobalIds[j], dofIndices[i], 1.0 );
        }//end for i
      }//end for j
    }//end for el
  }//end for nd

  mat->makeConsistent();
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


