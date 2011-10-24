#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/PIO.h"

#include "../MeshManager.h"
#include "../MeshVariable.h"

#include "vectors/Variable.h"

#include "libmesh_config.h"
#include "libmesh.h"
#include "mesh.h"
#include "cell_hex8.h"
#include "boundary_info.h"

#include "test_MeshGenerators.h"


void myTest(Parallel_Unit_Test *ut)
{

  std::string exeName("test_MakeSimpleMesh");
  std::string log_file = "log_" + exeName;

  AMP::AMP_MPI::initialize();
  AMP::PIO::logOnlyNodeZero(log_file);

  if(nodes == 1) {
    AMP::MeshManager::Adapter::shared_ptr meshAdapter = AMP::unit_test::ThreeElementLGenerator::getMesh();
    std::vector<std::vector<unsigned int> > elemNodeMap = AMP::unit_test::ThreeElementLGenerator::getElemNodeMap();
    std::vector<unsigned int> bndDofIndices = AMP::unit_test::ThreeElementLGenerator::getBndDofIndices();

    AMP::LinearAlgebra::Variable::shared_ptr var(new AMP::NodalScalarVariable("dummy"));
    AMP::DOFMap::shared_ptr dof_map = meshAdapter->getDOFMap(var);

    AMP::Mesh::MeshManager::Adapter::ElementIterator  el = meshAdapter->beginElement();
    AMP::Mesh::MeshManager::Adapter::ElementIterator  end_el = meshAdapter->endElement();
    AMP::Mesh::MeshManager::Adapter::ElementIterator  expectedEndEl = el;
    for(int i = 0; i < 3; i++) {
      expectedEndEl++; 
    }

    if(end_el != expectedEndEl) {
      ut.failure("end_el != expectedEndEl");
    }

    for(int i = 0; el != end_el; ++el, i++) {
      std::vector<unsigned int> dofIndices;
      dof_map->getDOFs(*el, dofIndices);
      unsigned int num_local_dofs = dofIndices.size();
      if(num_local_dofs != 8) {
        ut.numFails++;
        ut.failure("num_local_dofs != 8");
      }
      for(int j = 0; j < num_local_dofs; j++) {
        if(dofIndices[j] != elemNodeMap[i][j]) {
          ut.numFails++;
          ut.failure("dofIndices[j] != elemNodeMap[i][j]");
        }
      }
    }//end for el

    const short int boundaryId = 1;
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = meshAdapter->beginOwnedBoundary(boundaryId);
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = meshAdapter->endOwnedBoundary(boundaryId);
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator expectedEndBnd = bnd;
    for(int i = 0; i < 4; i++) {
      expectedEndBnd++; 
    }

    if(end_bnd != expectedEndBnd) {
      ut.failure("end_bnd != expectedEndBnd");
    }

    for(; bnd != end_bnd; ++bnd) {
      std::vector<unsigned int> bndGlobalIds;
      std::vector<unsigned int> singleton(1);
      singleton[0] = 0;
      dof_map->getDOFs(*bnd, bndGlobalIds, singleton);

      if(bndGlobalIds.size() != 1) {
        ut.numFails++;
        ut.failure("bndGlobalIds.size() != 1");
      }

      std::vector<unsigned int>::iterator pos = std::find(bndDofIndices.begin(), bndDofIndices.end(), bndGlobalIds[0]);

      if(pos == bndDofIndices.end()) {
        ut.numFails++;
        ut.failure("pos == bndDofIndices.end()");
      } else if(*pos != bndGlobalIds[0]) {
        ut.numFails++;
        ut.failure("*pos != bndGlobalIds[0]");
      }

      bndDofIndices.erase(pos);

    }//end for bnd

  }

  ut.passes(exeName);

  AMP::AMPManager::shutdown();

}


int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
        myTest(&ut);
    } catch (std::exception &err) {
        std::cout << "ERROR: While testing " <<argv[0] << err.what() << std::endl;
        ut.failure("Error while testing");
    } catch( ... ) {
        std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
        ut.failure("Error while testing");
    }

    ut.report ();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   



