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

#include "LibMeshAdapter.h"

#include "ReadTestMesh.h"

void RunLibMeshFaceStuff() {

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

  //const DofMap & dof_map = system.get_dof_map();

  FEType fe_type = system.variable_type(V_var);

  QGauss qrule (2, fe_type.default_quadrature_order());

  AutoPtr<FEBase> fe_2d  (FEBase::build(2, fe_type));
  AutoPtr<FEBase> fe_3d  (FEBase::build(3, fe_type));

  fe_2d->attach_quadrature_rule (&qrule);
  fe_3d->attach_quadrature_rule (&qrule);

  const std::vector<Point>& q_point2d = fe_2d->get_xyz();
  const std::vector<Point>& q_point3d = fe_3d->get_xyz();

  const std::vector<Real>& JxW2d = fe_2d->get_JxW();
  const std::vector<Real>& JxW3d = fe_3d->get_JxW();

  const std::vector<std::vector<Real> >& phi2d = fe_2d->get_phi();
  const std::vector<std::vector<Real> >& phi3d = fe_3d->get_phi();

  const std::vector<std::vector<RealGradient> >& dphi2d = fe_2d->get_dphi();
  const std::vector<std::vector<RealGradient> >& dphi3d = fe_3d->get_dphi();

  const std::vector<Point>& normal2d = (fe_2d->get_normals());
  const std::vector<Point>& normal3d = (fe_3d->get_normals());

  MeshBase::const_element_iterator       el     = mesh->local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh->local_elements_end(); 

  for( ; el != end_el; ++el) {
    const Elem* elem = *el;

    for(unsigned int s = 0; s < elem->n_sides(); s++) {
      AutoPtr<Elem> side (elem->build_side(s));
      fe_2d->reinit(side.get());
      fe_3d->reinit(elem, s);

      for(unsigned int ijk = 0; ijk < elem->n_nodes(); ijk++) {
        std::cout<<"Node # "<<ijk<<" = "<<elem->node(ijk)<<std::endl;
        //unsigned int pqr = elem->node(ijk);
        std::cout<<"Coordinates -- "<<(elem->point(ijk))(0)<<" "<<(elem->point(ijk))(1)<<" "<<(elem->point(ijk))(2)<<std::endl;
      }

      const short int bnd_id = (mesh->boundary_info)->boundary_id (elem, s);
      if(  bnd_id == 1 ){
        for(unsigned int qp = 0; qp < qrule.n_points(); qp++) 
        {
          for(int d = 0; d < 3; d++) {
            std::cout<<"Pt2d["<<qp<<"]["<<d<<"] = "<<q_point2d[qp](d)<<std::endl;
            std::cout<<"Pt3d["<<qp<<"]["<<d<<"] = "<<q_point3d[qp](d)<<std::endl;
          }

          std::cout<<"JxW2d["<<qp<<"] = "<<JxW2d[qp]<<std::endl;
          std::cout<<"JxW3d["<<qp<<"] = "<<JxW3d[qp]<<std::endl;

          std::cout<<"phi2d.size = "<<phi2d.size()<<std::endl;
          for(size_t l = 0; l < phi2d.size(); l++) {
            std::cout<<"phi2d["<<l<<"]["<<qp<<"] = "<<phi2d[l][qp]<<std::endl;
          }

          std::cout<<"phi3d.size = "<<phi3d.size()<<std::endl;
          for(size_t l = 0; l < phi3d.size(); l++) {
            std::cout<<"phi3d["<<l<<"]["<<qp<<"] = "<<phi3d[l][qp]<<std::endl;
          }

          std::cout<<"dphi2d.size = "<<dphi2d.size()<<std::endl;
          for(size_t l = 0; l < dphi2d.size(); l++) {
            for(int d = 0; d < 3; d++) {
              std::cout<<"dphi2d["<<l<<"]["<<qp<<"]["<<d<<"] = "<<dphi2d[l][qp](d)<<std::endl;
            }
          }

          std::cout<<"dphi3d.size = "<<dphi3d.size()<<std::endl;
          for(size_t l = 0; l < dphi3d.size(); l++) {
            for(int d = 0; d < 3; d++) {
              std::cout<<"dphi3d["<<l<<"]["<<qp<<"]["<<d<<"] = "<<dphi3d[l][qp](d)<<std::endl;
            }
          }

          std::cout<<"normal2d.size = "<<normal2d.size()<<std::endl;
          std::cout<<"normal3d.size = "<<normal3d.size()<<std::endl;
          std::cout<<"normal3d["<<qp<<"].Magnitude = "<<normal3d[qp].size()<<std::endl;
          for(int d = 0; d < 3; d++) {
            std::cout<<"normal3d["<<qp<<"]["<<d<<"] = "<<normal3d[qp](d)<<std::endl;
          }
          std::cout<<std::endl;
        }
      }
    }
  }
}


int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
        RunLibMeshFaceStuff();
        ut.passes("Ran test without crashing");
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

