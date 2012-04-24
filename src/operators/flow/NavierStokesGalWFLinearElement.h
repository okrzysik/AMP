
#ifndef included_AMP_NavierStokesGalWFLinearElement
#define included_AMP_NavierStokesGalWFLinearElement

#include <vector>

#include "boost/shared_ptr.hpp"

/* AMP files */
#include "operators/flow/FlowElement.h"
#include "operators/flow/NavierStokesConstants.h"

namespace AMP {
namespace Operator {

  class NavierStokesGalWFLinearElement : public FlowElement 
  {
    public :

      NavierStokesGalWFLinearElement(const boost::shared_ptr<ElementOperationParameters>& params)
        : FlowElement(params) {

          d_JxW = &(d_fe[0]->get_JxW());

          d_u_dphi = &(d_fe[0]->get_dphi());

          d_u_phi = &(d_fe[0]->get_phi());

          d_p_dphi = &(d_fe[1]->get_dphi());

          d_p_phi = &(d_fe[1]->get_phi());

          d_xyz = &(d_fe[0]->get_xyz());

        }

      ~NavierStokesGalWFLinearElement() {  }

      void setElementStiffnessMatrix( std::vector<std::vector<double> > & elementStiffnessMatrix )
      {
        d_elementStiffnessMatrix = &(elementStiffnessMatrix);
      }

      void setElementVectors( const std::vector<std::vector<double> > & elementInputVectors ) 
      {
        d_elementInputVectors = elementInputVectors;
      }

      void apply();


    protected :

      double d_density;

      double d_fmu;

      const std::vector<Real> *d_JxW; 

      const std::vector<std::vector<RealGradient> > *d_u_dphi; 

      const std::vector<std::vector<Real> > *d_u_phi; 

      const std::vector<std::vector<RealGradient> > *d_p_dphi; 

      const std::vector<std::vector<Real> > *d_p_phi; 

      const std::vector<Point> *d_xyz; 

      std::vector<std::vector<double> > d_elementInputVectors; 
      
      std::vector<std::vector<double> > *d_elementStiffnessMatrix;

    private :

  };

}
}

#endif


