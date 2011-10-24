
#include "NavierStokesGalWFLinearElement.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Operator {

  void NavierStokesGalWFLinearElement :: apply() 
  {
    const std::vector<Real> & JxW = (*d_JxW);

    const std::vector<std::vector<RealGradient> > & dphi = (*d_u_dphi);

    const std::vector<std::vector<Real> > & phi = (*d_u_phi);

    const std::vector<std::vector<Real> > & p_phi = (*d_p_phi);

    std::vector<std::vector<double> > & elementStiffnessMatrix = (*d_elementStiffnessMatrix);

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

      for(unsigned int k = 0; k < num_nodes; k++) {

        u += d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*phi[k][qp]; 
        v += d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*phi[k][qp]; 
        w += d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*phi[k][qp];

        dudx += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*dphi[k][qp](0));
        dudy += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*dphi[k][qp](1));
        dudz += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*dphi[k][qp](2));

        dvdx += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*dphi[k][qp](0));
        dvdy += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*dphi[k][qp](1));
        dvdz += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*dphi[k][qp](2));

        dwdx += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*dphi[k][qp](0));
        dwdy += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*dphi[k][qp](1));
        dwdz += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*dphi[k][qp](2));

      } //end for k

      d_density = d_transportModel->getDensity();
      d_fmu = d_transportModel->getViscosity();

      for(unsigned int i = 0; i < num_nodes; i++) {

          for(unsigned int j = 0; j < num_nodes; j++) {

            elementStiffnessMatrix[4*i    ][4*j     ] += JxW[qp]*( (dphi[i][qp]*dphi[j][qp])  +                                            // diffusion term
                                                                   (u*dphi[j][qp](0)+ v*dphi[j][qp](1)+w*dphi[j][qp](2))*phi[i][qp] +      // convection term
                                                                    dudx*phi[i][qp]*phi[j][qp]   );                                        // Newton term

            elementStiffnessMatrix[4*i    ][4*j + 1 ] += JxW[qp]*dudy*phi[i][qp]*phi[j][qp];                                               // Newton term
            
            elementStiffnessMatrix[4*i    ][4*j + 2 ] += JxW[qp]*dudz*phi[i][qp]*phi[j][qp];                                               // Newton term

            elementStiffnessMatrix[4*i + 1][4*j + 1 ] += JxW[qp]*( (dphi[i][qp]*dphi[j][qp]) +                                             // diffusion term
                                                                   (u*dphi[j][qp](0)+ v*dphi[j][qp](1)+w*dphi[j][qp](2))*phi[i][qp] +      // convection term
                                                                    dvdy*phi[i][qp]*phi[j][qp]);                                           // Newton term

            elementStiffnessMatrix[4*i + 1][4*j     ] += JxW[qp]*dvdx*phi[i][qp]*phi[j][qp];                                               // Newton term
            
            elementStiffnessMatrix[4*i + 1][4*j + 2 ] += JxW[qp]*dvdz*phi[i][qp]*phi[j][qp];                                               // Newton term

            elementStiffnessMatrix[4*i + 2][4*j + 2 ] += JxW[qp]*( (dphi[i][qp]*dphi[j][qp]) +                                             // diffusion term
                                                                   (u*dphi[j][qp](0)+ v*dphi[j][qp](1)+w*dphi[j][qp](2))*phi[i][qp] +      // convection term
                                                                    dwdz*phi[i][qp]*phi[j][qp]);                                           // Newton term

            elementStiffnessMatrix[4*i + 2][4*j     ] += JxW[qp]*dwdx*phi[i][qp]*phi[j][qp];                                               // Newton term
            
            elementStiffnessMatrix[4*i + 2][4*j + 1 ] += JxW[qp]*dwdy*phi[i][qp]*phi[j][qp];                                               // Newton term

          }//end for j

          for(unsigned int k = 0; k < d_elementInputVectors[NavierStokes::PRESSURE].size(); k++) {
          
            elementStiffnessMatrix[4*i    ][4*k + 3] += JxW[qp]*(p_phi[k ][qp]*dphi[i][qp](0));
            
            elementStiffnessMatrix[4*i + 1][4*k + 3] += JxW[qp]*(p_phi[k ][qp]*dphi[i][qp](1));
            
            elementStiffnessMatrix[4*i + 2][4*k + 3] += JxW[qp]*(p_phi[k ][qp]*dphi[i][qp](2));
          
          }//end for k
      
      }//end for i

    }//end for qp

  }

}
}



