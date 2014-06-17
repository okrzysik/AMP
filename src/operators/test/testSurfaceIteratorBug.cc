#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "mpi.h"

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "ampmesh/Mesh.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorSelector.h"

#include "libmesh/fe_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/elem.h"
#include "libmesh/quadrature.h"
#include "libmesh/face_quad4.h"

#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/string_to_enum.h"


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

  // Get the Mesh database and create the mesh parameters
  boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
  boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
  params->setComm(globalComm);

  // Create the meshes from the input database
  boost::shared_ptr<AMP::Mesh::Mesh> mesh = AMP::Mesh::Mesh::buildMesh(params);

  // Create a nodal scalar vector
  AMP::LinearAlgebra::Variable::shared_ptr var(new  AMP::LinearAlgebra::Variable("myVar")); 
  AMP::Discretization::DOFManager::shared_ptr nodalScalarDOF = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,1,true);  
  AMP::LinearAlgebra::Vector::shared_ptr vec = AMP::LinearAlgebra::createVector( nodalScalarDOF, var, true );
  vec->zero();

  libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>("FIRST");
  libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>("LAGRANGE");

  AMP::Mesh::MeshIterator bnd = mesh->getBoundaryIDIterator( AMP::Mesh::Face, surfaceId, 0 );
  AMP::Mesh::MeshIterator end_bnd = bnd.end();

  bool volume_passes = true;
  std::vector<size_t> dofs;
  while ( bnd!=end_bnd ) {

    std::vector<AMP::Mesh::MeshElement> nodes = bnd->getElements(AMP::Mesh::Vertex);
    std::vector<size_t> bndGlobalIds;
    for (size_t i=0; i<nodes.size(); i++) {
        nodalScalarDOF->getDOFs(nodes[i].globalID(), dofs);
        for (size_t j=0; j<dofs.size(); j++)
            bndGlobalIds.push_back(dofs[j]);
    }
    
    // Some basic checks
    assert(bndGlobalIds.size() == 4);
    //assert((bnd->getElem()).default_order() == feTypeOrder);
    assert(bnd->elementType()==AMP::Mesh::Face);

    // Create the libmesh element
    // Note: This must be done inside the loop because libmesh's reinit function doesn't seem to work properly
    boost::shared_ptr < ::FEType > feType( new ::FEType(feTypeOrder, feFamily) );
    boost::shared_ptr < ::FEBase > fe( (::FEBase::build(2, (*feType))).release() );
    const std::vector<std::vector<Real> > &phi = fe->get_phi();
    const std::vector<Real> &djxw = fe->get_JxW();
    libMeshEnums::QuadratureType qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>("QGAUSS");
    libMeshEnums::Order qruleOrder = feType->default_quadrature_order();
    boost::shared_ptr < ::QBase > qrule( (::QBase::build(qruleType, 2, qruleOrder)).release() );
    fe->attach_quadrature_rule( qrule.get() );
    ::Elem* currElemPtr = new ::Quad4;
    for(size_t i=0; i<nodes.size(); i++) {
        std::vector<double> pt = nodes[i].coord();
        currElemPtr->set_node(i) = new ::Node(pt[0], pt[1], pt[2], i);
    }
    fe->reinit( currElemPtr );

    // Check the volume
    double vol1 = 0.0;
    for(unsigned int qp = 0; qp < qrule->n_points(); qp++)
        vol1 += djxw[qp];
    double vol2 = bnd->volume();
    if( fabs(vol1-vol2) > (1.0e-8*vol2) ) {
        std::cout<<"Volume 1 = "<<std::setprecision(15)<<vol1<<std::endl;
        std::cout<<"Volume 2 = "<<std::setprecision(15)<<vol2<<std::endl<<std::endl;
        volume_passes = false;
    }

    // Fill the surface vector
    if(setConstantValue) {
       std::vector<double> vals(bndGlobalIds.size(), 100.0);
       vec->addValuesByGlobalID(bndGlobalIds.size(), &(bndGlobalIds[0]), &(vals[0]));
    } else {
       std::vector<double> vals(bndGlobalIds.size(), 0.0);
       for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
         for(unsigned int qp = 0; qp < qrule->n_points(); qp++) {
           assert(djxw[qp] >= 0.0);
           assert(phi[i][qp] >= 0.0);
           vals[i] += (djxw[qp]*phi[i][qp]*100.0);
         }//end qp
       }//end i
       vec->addValuesByGlobalID(bndGlobalIds.size(), &(bndGlobalIds[0]), &(vals[0]));
    }

    // Destroy the libmesh element
    for(size_t i=0; i<nodes.size(); i++) {
        delete (currElemPtr->get_node(i));
        currElemPtr->set_node(i) = NULL;
    }//end for j
    delete currElemPtr;
    currElemPtr = NULL;

    ++bnd;
  }//end for bnd

  if ( volume_passes == true )
    ut->passes("Volume passes");
  else
    ut->failure("Volume passes");

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );

  double l2Norm = vec->L2Norm();
  std::cout<<"size = "<<vec->getGlobalSize()<<std::endl;
  std::cout<<"L2 Norm = "<<std::setprecision(15)<<l2Norm<<std::endl;

  if ( AMP::Utilities::approx_equal(l2Norm,0.00106829009941852) )
     ut->passes("L2 Norm has expected value");
  else
     ut->failure("L2 Norm has expected value");
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testSurfaceIteratorBug";

  try {
    myTest(&ut, exeName);
    ut.passes(exeName);
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



