#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include <iostream>
#include <cstdlib>

#include "boost/shared_ptr.hpp"

#include "mesh.h"
#include "mesh_communication.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "linear_implicit_system.h"
#include "elem.h"
#include "boundary_info.h"


#include "utils/ReadTestMesh.h"

/* Libmesh files */
#include "fe_type.h"
#include "fe_base.h"
#include "elem.h"
#include "quadrature.h"

#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"
#include "string_to_enum.h"

//Using mesh and function calls from testLibmeshFaceStuff.cc

void calculateGrad(AMP::UnitTest *ut)
{
  const unsigned int mesh_dim = 3;
  boost::shared_ptr< Mesh > mesh(new Mesh(mesh_dim));

  AMP::readTestMesh("distortedElementMesh", mesh);

  MeshCommunication().broadcast(*(mesh.get()));

  mesh->prepare_for_use(false);

  EquationSystems equation_systems (*(mesh.get()));

  LinearImplicitSystem & system = 
    equation_systems.add_system<LinearImplicitSystem> ("Poisson");

  system.add_variable ("V", FIRST);
  equation_systems.init ();

  const unsigned int V_var = system.variable_number ("V");

  FEType fe_type = system.variable_type(V_var);

  QGauss qrule (3, fe_type.default_quadrature_order());

  AutoPtr<FEBase> fe_3d  (FEBase::build(3, fe_type));
  fe_3d->attach_quadrature_rule (&qrule);
  //const std::vector<Point>& q_point3d = fe_3d->get_xyz();
  
  const std::vector<std::vector<Real> >& dphi3d = fe_3d->get_dphidx();

  MeshBase::const_element_iterator       el     = mesh->local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->local_elements_end(); 
  std::cout << "Entering Element Iyerator" << std::endl;
  for( ; el != end_el; ++el) {
    const Elem* elem = *el;

    fe_3d->reinit(elem);
    std::vector<Point> coordinates = fe_3d->get_xyz();
    std::vector<double>  computedAtGauss(qrule.n_points(), 0.0);
    std::cout << "Entering Gauss Point loop : "<< qrule.n_points()<< std::endl;
    for(unsigned int qp = 0; qp < qrule.n_points(); qp++) 
    {
      std::cout<<"dphidx.size = "<<dphi3d.size()<<std::endl;
      for(size_t l = 0; l < dphi3d.size(); l++) {
          std::cout<<"dphidx["<<l<<"]["<<qp<<"]  = "<<dphi3d[l][qp]<<std::endl;
      }
    }
  }
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::cout << "Entering main" << std::endl;
    try
    {
      calculateGrad(&ut);
    } catch (std::exception &err) {
      std::cout << "ERROR:While testing " << argv[0]<<err.what()<<std::endl;
      ut.failure("ERROR: While testing");
    } catch(...){
      std::cout   << "ERROR: While testing " << argv[0] << "An unknown exception was thrown." << std::endl;
      ut.failure("ERROR: While testing");
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
