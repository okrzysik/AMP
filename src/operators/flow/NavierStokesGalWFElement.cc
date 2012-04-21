
#include "NavierStokesGalWFElement.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Operator {

  void NavierStokesGalWFElement :: initTransportModel() 
  {
    (d_fe[0])->reinit(d_elem);
    (d_fe[1])->reinit(d_elem);

  }

  void NavierStokesGalWFElement :: apply() 
  {
    const std::vector<Real> & JxW = (*d_JxW);

    const std::vector<std::vector<RealGradient> > & u_dphi = (*d_u_dphi);

    const std::vector<std::vector<Real> > & u_phi = (*d_u_phi);

    const std::vector<std::vector<RealGradient> > & p_dphi = (*d_p_dphi);

    const std::vector<std::vector<Real> > & p_phi = (*d_p_phi);

    std::vector<std::vector<double> > & elementInputVectors = d_elementInputVectors;

    std::vector<double> & elementOutputVector = (*d_elementOutputVector);

    (d_fe[0])->reinit(d_elem);
    (d_fe[1])->reinit(d_elem);


    const unsigned int num_nodes = d_elem->n_nodes();

    for(unsigned int qp = 0; qp < (d_qrule[0])->n_points(); qp++) {

      double u = 0;
      double v = 0;
      double w = 0;

      double dudx = 0;
      double dudy = 0;
      double dudz = 0;
      double dvdx = 0;
      double dvdy = 0;
      double dvdz = 0;
      double dwdx = 0;
      double dwdy = 0;
      double dwdz = 0;

      double p = 0;
      double dpdx = 0;
      double dpdy = 0;
      double dpdz = 0;

      for(unsigned int k = 0; k < num_nodes; k++) {

        u += elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*u_phi[k][qp]; 
        v += elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*u_phi[k][qp]; 
        w += elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*u_phi[k][qp];
        p += elementInputVectors[NavierStokes::PRESSURE][k]*p_phi[k][qp];


        dudx += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*u_dphi[k][qp](0));
        dudy += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*u_dphi[k][qp](1));
        dudz += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*u_dphi[k][qp](2));

        dvdx += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*u_dphi[k][qp](0));
        dvdy += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*u_dphi[k][qp](1));
        dvdz += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*u_dphi[k][qp](2));

        dwdx += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*u_dphi[k][qp](0));
        dwdy += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*u_dphi[k][qp](1));
        dwdz += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*u_dphi[k][qp](2));

        dpdx += (elementInputVectors[NavierStokes::PRESSURE][k ]*p_dphi[k][qp](0));
        dpdy += (elementInputVectors[NavierStokes::PRESSURE][k ]*p_dphi[k][qp](1));
        dpdz += (elementInputVectors[NavierStokes::PRESSURE][k ]*p_dphi[k][qp](2));

      } //end for k

      d_density = d_transportModel->getDensity();
      d_fmu = d_transportModel->getViscosity();

      for(unsigned int j = 0; j < num_nodes; j++) {

        elementOutputVector[(4*j)     ] += JxW[qp]*( d_alpha_conv * ( u*dudx + v*dudy + w*dudz ) * u_phi[j][qp] + d_density * dpdx * p_phi[j][qp] + 
           d_alpha_diff * d_fmu *( dudx * u_dphi[j][qp](0) + dudy * u_dphi[j][qp](1) + dudz * u_dphi[j][qp](2) ) );

        elementOutputVector[(4*j) + 1 ] += JxW[qp]*( d_alpha_conv * ( u*dvdx + v*dvdy + w*dvdz ) * u_phi[j][qp] + d_density * dpdy * p_phi[j][qp] + 
           d_alpha_diff * d_fmu *( dvdx * u_dphi[j][qp](0) + dvdy * u_dphi[j][qp](1) + dvdz * u_dphi[j][qp](2) ) );

        elementOutputVector[(4*j) + 2 ] += JxW[qp]*( d_alpha_conv * ( u*dwdx + v*dwdy + w*dwdz ) * u_phi[j][qp] + d_density * dpdz * p_phi[j][qp] + 
           d_alpha_diff * d_fmu *( dwdx * u_dphi[j][qp](0) + dwdy * u_dphi[j][qp](1) + dwdz * u_dphi[j][qp](2) ) );

        elementOutputVector[(4*j) + 3 ] += JxW[qp]*( ( dudx + dudy + dudz ) * u_phi[j][qp] );

      }

    }//end for qp

  }

}
}

