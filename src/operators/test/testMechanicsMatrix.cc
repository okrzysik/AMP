
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

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"
#include "matrices/MatrixBuilder.h"

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

  AMP::Discretization::DOFManager::shared_ptr nodal3VectorDOF = AMP::Discretization::simpleDOFManager::create(meshAdapter,AMP::Mesh::Vertex, 1, 3,true);  
  boost::shared_ptr<AMP::LinearAlgebra::Variable> inpVariable (new AMP::LinearAlgebra::Variable("InputVar") );
  boost::shared_ptr<AMP::LinearAlgebra::Variable> outVariable (new AMP::LinearAlgebra::Variable("OutVar") );

  AMP::LinearAlgebra::Vector::shared_ptr inVec  = AMP::LinearAlgebra::createVector(nodal3VectorDOF , inpVariable, true);
  AMP::LinearAlgebra::Vector::shared_ptr outVec = AMP::LinearAlgebra::createVector(nodal3VectorDOF , outVariable, true);

  AMP::LinearAlgebra::Matrix::shared_ptr mat = AMP::LinearAlgebra::createMatrix(inVec, outVec);

  mat->setScalar(1);

  //------------------------------
  //DirichletCorrection-----------
  //------------------------------
  
  AMP::Mesh::MeshIterator bnd = meshAdapter->getIDsetIterator( AMP::Mesh::Vertex, 2, 0 );
  AMP::Mesh::MeshIterator end_bnd = bnd.end();

  for( ; bnd != end_bnd; ++bnd) {
    std::vector<size_t> bndDofIds;
    nodal3VectorDOF->getDOFs(bnd->globalID(), bndDofIds);

    std::vector< AMP::Mesh::MeshElement::shared_ptr > neighbors = bnd->getNeighbors();
    for(unsigned int i = 0; i < neighbors.size(); ++i) {
      AMP_ASSERT((*(neighbors[i])) != (*bnd));
    }//end for i

    for(unsigned int j = 0; j < 3; ++j) {
      for(unsigned int i = 0; i < bndDofIds.size(); ++i) {
        if(j == i) {
          mat->setValueByGlobalID ( bndDofIds[i], bndDofIds[i], 1.0 );
        } else {
          mat->setValueByGlobalID ( bndDofIds[j], bndDofIds[i], 0.0 );
          mat->setValueByGlobalID ( bndDofIds[i], bndDofIds[j], 0.0 );
        }
      }//end for i
      for(size_t n = 0; n < neighbors.size(); ++n) {
        std::vector<size_t> nhDofIds;
        nodal3VectorDOF->getDOFs(neighbors[n]->globalID(), nhDofIds);
        for(unsigned int i = 0; i < nhDofIds.size(); ++i) {
          mat->setValueByGlobalID ( bndDofIds[j], nhDofIds[i], 0.0 );
          mat->setValueByGlobalID ( nhDofIds[i], bndDofIds[j], 0.0 );
        }//end for i
      }//end for n
    }//end for j
  }//end for bnd

  size_t locSize = inVec->getLocalSize();
  size_t globSize = inVec->getGlobalSize();
  size_t locStartId = inVec->getLocalStartID();

  AMP::plog<<" locStartID = "<<locStartId<<" locSize = "<<locSize<<" globSize = "<<globSize<<std::endl;

  AMP::Mesh::MeshIterator nd = meshAdapter->getIterator(AMP::Mesh::Vertex, 0);
  AMP::Mesh::MeshIterator end_nd = nd.end();
  for(int cnt = 0; nd != end_nd; ++nd, ++cnt) {
    std::vector<size_t> ndDofIds;
    nodal3VectorDOF->getDOFs(nd->globalID(), ndDofIds);
    std::vector<double> pt = nd->coord();
    AMP::plog<<std::endl<<" locNdCnt = "<<cnt<<" Pt: "<<(pt[0])<<" : "<<(pt[1])<<" : "<<(pt[2])<<std::endl;
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



