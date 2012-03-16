
#include <iostream>
#include <string>
#include <cassert>
#include <fstream>

#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

/* AMP files */
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

void myTest(AMP::UnitTest *ut )
{
  std::string exeName("testMechanicsMatrix");
  std::string input_file = "input_" + exeName;

  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  int npes = globalComm.getSize();

  char outFile[256];
  sprintf(outFile, "outMechMat_%d", npes);

  AMP::PIO::logAllNodes(outFile);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( mgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ( "bar" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator =
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "MechanicsBVPOperator", input_db, elementPhysicsModel));

  AMP::LinearAlgebra::Matrix::shared_ptr mat = bvpOperator->getMatrix();
  AMP::LinearAlgebra::Vector::shared_ptr vec = mat->getLeftVector();
  AMP::LinearAlgebra::Variable::shared_ptr var = bvpOperator->getOutputVariable();

  AMP::Mesh::DOFMap::shared_ptr dofMap = meshAdapter->getDOFMap(var);

  size_t locSize = vec->getLocalSize();
  size_t globSize = vec->getGlobalSize();
  size_t locStartId = vec->getLocalStartID();

  AMP::plog<<" locStartID = "<<locStartId<<" locSize = "<<locSize<<" globSize = "<<globSize<<std::endl;

  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator nd = meshAdapter->beginOwnedNode();
  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator end_nd = meshAdapter->endOwnedNode();
  for(int cnt = 0; nd != end_nd; ++nd, ++cnt) {
    std::vector<unsigned int> ndDofIds;
    std::vector<unsigned int> empty;
    dofMap->getDOFs(*nd, ndDofIds, empty);
    AMP::plog<<std::endl<<" locNdCnt = "<<cnt<<" Pt: "<<(nd->x())<<" : "<<(nd->y())<<" : "<<(nd->z())<<std::endl;
    for(int i = 0; i < ndDofIds.size(); ++i) {
      std::vector<unsigned int> cols;
      std::vector<double> vals;
      mat->getRowByGlobalID(ndDofIds[i], cols, vals);
      AMP::plog<<std::endl<<" row = "<<(ndDofIds[i])<<" NumCols = "<<(cols.size())<<std::endl;
      for(int j = 0; j < cols.size(); ++j) {
        AMP::plog<<" col("<<j<<") = "<<(cols[j])<<" Val = "<<(vals[j])<<std::endl;
      }//end j
    }//end i
  }//end nd

  ut->passes(exeName);

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



