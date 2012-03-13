
#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

//#include "ContactSearchUtils.h"

void myTest(AMP::UnitTest *ut, std::string exeName) {
  
#if 0
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  const double precision = 1.0e-12;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  const std::string masterMeshName = input_db->getString("MasterMesh");
  const std::string slaveMeshName = input_db->getString("SlaveMesh");
  const unsigned int slaveId = input_db->getInteger("SlaveId");
  const unsigned int masterId = input_db->getInteger("MasterId");
  const unsigned int rgDim = input_db->getInteger("Num1Dcells");


  AMP::Mesh::MeshManagerParameters::shared_ptr meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter = manager->getMesh ( masterMeshName );
  AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter = manager->getMesh ( slaveMeshName );

  double minXYZ[3];
  double maxXYZ[3];

  computeRGboundingBox(precision, masterMeshAdapter, minXYZ, maxXYZ);

  double rgH[3];
  std::vector<std::vector<size_t> > rg2ElemMap;

  computeRG2ElemMap(precision, rgDim, masterMeshAdapter, minXYZ, maxXYZ, rg2ElemMap, rgH);

  std::vector<size_t> slaveNodes;
  std::vector<size_t> slave2MasterElem;

  computeSlave2MasterElem(slaveId, masterId, precision, rgH, rgDim, 
      slaveMeshAdapter, masterMeshAdapter, minXYZ, maxXYZ, rg2ElemMap, slaveNodes, slave2MasterElem);

  std::vector<std::vector<size_t> > slave2MasterNodes;
  std::vector<std::vector<double> > slave2MasterFactors;

  computeSlave2MasterNodes(precision, slaveId, masterId, slaveMeshAdapter, masterMeshAdapter,
      slaveNodes, slave2MasterElem, slave2MasterNodes, slave2MasterFactors);

#endif

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testContactSearch";

  AMP_ERROR("Not yet converted!");

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


