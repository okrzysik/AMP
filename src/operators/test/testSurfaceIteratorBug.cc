
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

#include "fe_type.h"
#include "fe_base.h"
#include "elem.h"
#include "quadrature.h"

#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"
#include "string_to_enum.h"

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  int surfaceId = input_db->getInteger("SurfaceId");
  bool setConstantValue = input_db->getBool("SetConstantValue");

  AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( mgrParams ) );
  AMP::Mesh::MeshAdapter::shared_ptr mesh = manager->getMesh("pellet");

  AMP::LinearAlgebra::Variable::shared_ptr var(new AMP::Mesh::NodalScalarVariable("myVar", mesh)); 
  AMP::Mesh::DOFMap::shared_ptr dof_map = mesh->getDOFMap(var);

  AMP::LinearAlgebra::Vector::shared_ptr vec = mesh->createVector(var);
  vec->zero();

  libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>("FIRST");
  libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>("LAGRANGE");
  boost::shared_ptr < ::FEType > feType( new ::FEType(feTypeOrder, feFamily) );

  libMeshEnums::QuadratureType qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>("QGAUSS");
  libMeshEnums::Order qruleOrder = feType->default_quadrature_order();
  boost::shared_ptr < ::QBase > qrule( (::QBase::build(qruleType, 2, qruleOrder)).release() );

  boost::shared_ptr < ::FEBase > fe( (::FEBase::build(2, (*feType))).release() );
  fe->attach_quadrature_rule( qrule.get() );

  const std::vector<std::vector<Real> > &phi = fe->get_phi();
  const std::vector<Real> &djxw = fe->get_JxW();

  AMP::Mesh::MeshManager::Adapter::BoundarySideIterator bnd = mesh->beginSideBoundary( surfaceId );
  AMP::Mesh::MeshManager::Adapter::BoundarySideIterator end_bnd = mesh->endSideBoundary( surfaceId );

  for( ; bnd != end_bnd; ++bnd) {
    std::vector<unsigned int> bndGlobalIds;
    dof_map->getDOFs(*bnd, bndGlobalIds, 0);

    assert(bndGlobalIds.size() == 4);

    if(setConstantValue) {
      std::vector<double> vals(bndGlobalIds.size(), 100.0);
      vec->addValuesByGlobalID(bndGlobalIds.size(), (int*)(&(bndGlobalIds[0])), &(vals[0]));
    } else {
      fe->reinit ( &(bnd->getElem()) );
      std::vector<double> vals(bndGlobalIds.size(), 0.0);
      for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
        for(unsigned int qp = 0; qp < qrule->n_points(); qp++) {
          vals[i] +=  (djxw[qp]*phi[i][qp]*100.0);
        }//end qp
      }//end i
      vec->addValuesByGlobalID(bndGlobalIds.size(), (int*)(&(bndGlobalIds[0])), &(vals[0]));
    }
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



