
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

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator =
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "MechanicsBVPOperator", input_db, elementPhysicsModel));

  AMP::LinearAlgebra::Matrix::shared_ptr mat = bvpOperator->getMatrix();
  AMP::LinearAlgebra::Vector::shared_ptr vec = mat->getLeftVector();

  AMP::Discretization::DOFManager::shared_ptr dofMap = vec->getDOFManager();

  size_t locSize = vec->getLocalSize();
  size_t globSize = vec->getGlobalSize();
  size_t locStartId = vec->getLocalStartID();

  AMP::plog<<" locStartID = "<<locStartId<<" locSize = "<<locSize<<" globSize = "<<globSize<<std::endl;

  AMP::Mesh::MeshIterator nd = meshAdapter->getIterator(AMP::Mesh::Vertex, 0);
  AMP::Mesh::MeshIterator end_nd = nd.end();
  for(int cnt = 0; nd != end_nd; ++nd, ++cnt) {
    std::vector<size_t> ndDofIds;
    dofMap->getDOFs(nd->globalID(), ndDofIds);
    std::vector<double> pt = nd->coord();
    AMP::plog<<" locNdCnt = "<<cnt<<" Pt: "<<(pt[0])<<" : "<<(pt[1])<<" : "<<(pt[2])<<std::endl;
    for(int i = 0; i < ndDofIds.size(); ++i) {
      std::vector<unsigned int> cols;
      std::vector<double> vals;
      mat->getRowByGlobalID(ndDofIds[i], cols, vals);
      AMP::plog<<" row = "<<(ndDofIds[i])<<" NumCols = "<<(cols.size())<<std::endl;
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



