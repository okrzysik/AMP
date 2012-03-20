
#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>
#include <cstdlib>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/AMP_MPI.h"

#include "vectors/PetscVector.h"

#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshAdapter.h"
#include "ampmesh/MeshVariable.h"

#include "libmesh.h"

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
  AMP::Mesh::MeshAdapter::shared_ptr mesh = manager->getMesh("TestMesh");

  const unsigned int DOFsPerNode = 1;
  AMP::LinearAlgebra::Variable::shared_ptr var (new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, DOFsPerNode>("var", mesh) ) ; 
  AMP::Mesh::DOFMap::shared_ptr dof_map = mesh->getDOFMap(var);

  boost::shared_ptr<AMP::LinearAlgebra::Matrix> mat =  mesh->createMatrix ( var, var );
  mat->setScalar(2);
  AMP::LinearAlgebra::Vector::shared_ptr vec = mesh->createVector(var);
  vec->zero();


  AMP::Mesh::MeshManager::Adapter::BoundaryNodeIterator bnd = mesh->beginBoundary( 2 );
  AMP::Mesh::MeshManager::Adapter::BoundaryNodeIterator end_bnd = mesh->endBoundary( 2 );
  std::vector<size_t> bnd_dofs;
  for( ; bnd != end_bnd; ++bnd) {
    AMP::Mesh::MeshManager::Adapter::NodeElementIterator el =  mesh->beginElementForNode ( *bnd );
    AMP::Mesh::MeshManager::Adapter::NodeElementIterator end_el = mesh->endElementForNode ( *bnd );

    std::vector<unsigned int> bndGlobalIds;
    std::vector<unsigned int> d_dofIds;
    dof_map->getDOFs(*bnd, bndGlobalIds, d_dofIds);
    for (size_t i=0; i<bndGlobalIds.size(); i++)
      bnd_dofs.push_back(bndGlobalIds[i]);

    for( ; el != end_el; ++el) {
      std::vector<unsigned int> dofIndices;
      dof_map->getDOFs(*el, dofIndices);

      for(unsigned int j = 0; j < bndGlobalIds.size(); j++) {
        for(unsigned int i = 0; i < dofIndices.size(); i++) {
          if(bndGlobalIds[j] == dofIndices[i]) {
            mat->setValueByGlobalID ( bndGlobalIds[j], bndGlobalIds[j], 1.0 );
          } else {
            mat->setValueByGlobalID ( bndGlobalIds[j], dofIndices[i] , 0.0 );
            mat->setValueByGlobalID ( dofIndices[i], bndGlobalIds[j] , 0.0 );
          }
        }//end for i
      }//end for j
    }//end for el
  }//end for bnd

  mat->makeConsistent();
  //-------------------------------------------

  AMP::Utilities::quicksort(bnd_dofs);
  AMP::plog << "Boundry DOFS:" << std::endl;
  for (size_t i=0; i<bnd_dofs.size(); i++)
     AMP::plog << "   " << bnd_dofs[i] << std::endl;

  size_t locSize = vec->getLocalSize();
  size_t globSize = vec->getGlobalSize();
  size_t locStartId = vec->getLocalStartID();

  AMP::plog<<" locStartID = "<<locStartId<<" locSize = "<<locSize<<" globSize = "<<globSize<<std::endl;

  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator  nd = mesh->beginOwnedNode();
  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator  end_nd =  mesh->endOwnedNode();

  for(int cnt = 0; nd != end_nd; ++nd, ++cnt) {

    std::vector<unsigned int> ndDofIds;
    std::vector<unsigned int> d_dofIds;
    d_dofIds.resize(0);
    dof_map->getDOFs(*nd, ndDofIds, d_dofIds);

    AMP::plog<<std::endl<<" locNdCnt = "<<cnt<<" Pt: "<<((*nd).x())<<" : "<<((*nd).y())<<" : "<<((*nd).z())<<std::endl;

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



