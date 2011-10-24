
#include "ConsMassGalWFLinearElement.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Operator {

  void ConsMassGalWFLinearElement :: apply() 
  {
    const std::vector<Real> & JxW = (*d_JxW);

    const std::vector<std::vector<RealGradient> > & u_dphi = (*d_u_dphi);

    const std::vector<std::vector<Real> > & p_phi = (*d_p_phi);

    std::vector<std::vector<double> > & elementStiffnessMatrix = (*d_elementStiffnessMatrix);

    (d_fe[0])->reinit(d_elem);
    (d_fe[1])->reinit(d_elem);

    for(unsigned int qp = 0; qp < (d_qrule[0])->n_points(); qp++) {

      for(unsigned int i = 0; i < u_dphi.size(); i++) {
        for(unsigned int dr = 0; dr < 3; dr++) {
          for(unsigned int j = 0; j < p_phi.size(); j++) {

            elementStiffnessMatrix[3*i + dr][j] += JxW[qp]*( (p_phi[j][qp]*u_dphi[i][qp](dr)) );                                 

          }//end for j
        }
      }//end for i

    }//end for qp

  }

}
}



