
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "mpi.h"

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshAdapter.h"
#include "ampmesh/MeshVariable.h"

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  int surfaceId = input_db->getInteger("SurfaceId");

  AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( mgrParams ) );
  AMP::Mesh::MeshAdapter::shared_ptr mesh = manager->getMesh("pellet");

  AMP::LinearAlgebra::Variable::shared_ptr var(new AMP::Mesh::NodalScalarVariable("myVar", mesh)); 
  AMP::Mesh::DOFMap::shared_ptr dof_map = mesh->getDOFMap(var);

  AMP::LinearAlgebra::Vector::shared_ptr vec = mesh->createVector(var);
  vec->zero();

  AMP::Mesh::MeshManager::Adapter::BoundarySideIterator bnd = mesh->beginSideBoundary( surfaceId );
  AMP::Mesh::MeshManager::Adapter::BoundarySideIterator end_bnd = mesh->endSideBoundary( surfaceId );

  for( ; bnd != end_bnd; ++bnd) {
    std::vector<unsigned int> bndGlobalIds;
    dof_map->getDOFs(*bnd, bndGlobalIds, 0);

    assert(bndGlobalIds.size() == 4);

    std::vector<double> vals(bndGlobalIds.size(), 10.0);
    vec->addValuesByGlobalID(bndGlobalIds.size(), (int*)(&(bndGlobalIds[0])), &(vals[0]));
  }//end for bnd

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );

  double l2Norm = vec->L2Norm();
  std::cout<<"L2 Norm = "<<std::setprecision(15)<<l2Norm<<std::endl;

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testSurfaceIteratorBug";

  try {
    myTest(&ut, exeName);
  } catch (std::exception &err) {
    std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
    ut.failure("ERROR: While testing");
  } catch( ... ) {
    std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
    ut.failure("ERROR: While testing");
  }

  //  ut.report();

  int num_failed = ut.NumFailGlobal();
  AMP::AMPManager::shutdown();
  return num_failed;
}  



