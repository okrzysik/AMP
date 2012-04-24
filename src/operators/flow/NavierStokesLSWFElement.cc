
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

    d_fmu =  d_transportModel->getViscosity();
    
    double d_Re  =  d_transportModel->getReynoldsNumber();

    const unsigned int num_nodes = d_elem->n_nodes();

    for(unsigned int qp = 0; qp < (d_qrule[0])->n_points(); qp++) {

      double u = 0;
      double v = 0;
      double w = 0;

      double txx = 0;
      double tyy = 0;
      double tzz = 0;
      double txy = 0;
      double tyz = 0;
      double txz = 0;

      double dtxxdx = 0;
      double dtyydy = 0;
      double dtzzdz = 0;

      double dtxydx = 0;
      double dtxydy = 0;

      double dtyzdy = 0;
      double dtyzdz = 0;

      double dtxzdx = 0;
      double dtxzdz = 0;

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

        txx += elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) + 0]*u_phi[k][qp]; 
        tyy += elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) + 1]*u_phi[k][qp]; 
        tzz += elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) + 2]*u_phi[k][qp];
        
        dtxxdx += (elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) + 0]*u_dphi[k][qp](0));
        dtyydy += (elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) + 1]*u_dphi[k][qp](1));
        dtzzdz += (elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) + 2]*u_dphi[k][qp](2));

        txy += elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) + 0]*u_phi[k][qp]; 
        tyz += elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) + 1]*u_phi[k][qp]; 
        txz += elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) + 2]*u_phi[k][qp];
        
        dtxydx += (elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) + 0]*u_dphi[k][qp](0));
        dtyzdy += (elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) + 1]*u_dphi[k][qp](1));
        dtxzdz += (elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) + 2]*u_dphi[k][qp](2));

        dtxydy += (elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) + 0]*u_dphi[k][qp](1));
        dtyzdz += (elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) + 1]*u_dphi[k][qp](2));
        dtxzdx += (elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) + 2]*u_dphi[k][qp](0));
      } //end for k
      
      d_density = d_transportModel->getDensity();
      d_fmu = d_transportModel->getViscosity();

      std::vector<double> error(10);
      error[0] = (dudx + dvdy + dwdz ) ;
      error[1] = (u * dudx + v * dudy + w * dudz ) + dpdx + ( dtxxdx + dtxydy + dtxzdz) ;
      error[2] = (u * dvdx + v * dvdy + w * dvdz ) + dpdy + ( dtxydx + dtyydy + dtyzdz) ;
      error[3] = (u * dwdx + v * dwdy + w * dwdz ) + dpdz + ( dtxzdx + dtyzdy + dtzzdz) ;
      error[4] = txx + 2.0 * (d_fmu/d_Re) * dudx; 
      error[5] = tyy + 2.0 * (d_fmu/d_Re) * dvdy; 
      error[6] = tzz + 2.0 * (d_fmu/d_Re) * dwdz; 
      error[7] = txy + (d_fmu/d_Re) * ( dudy + dvdx ); 
      error[8] = tyz + (d_fmu/d_Re) * ( dvdz + dwdy ); 
      error[9] = txz + (d_fmu/d_Re) * ( dwdx + dudz );  

      for(unsigned int j = 0; j < num_nodes; j++) {
        for(unsigned int i = 0; i < error.size() ; j++) {
          elementOutputVector[(4*j) + i ] += JxW[qp] * ( error[i] * error[i] ) * u_phi[j][qp];
        }
      }

    }//end for qp

  }

}
}

