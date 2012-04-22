
#ifndef included_AMP_NavierStokesLSWFElement
#define included_AMP_NavierStokesLSWFElement

#include <vector>

#include "boost/shared_ptr.hpp"

/* AMP files */
#include "operators/flow/FlowElement.h"
#include "operators/flow/NavierStokesConstants.h"

namespace AMP {
namespace Operator {

  class NavierStokesLSWFElement : public FlowElement 
  {
    public :

      /**
        Constructor.
        */
      NavierStokesLSWFElement(const boost::shared_ptr<ElementOperationParameters>& params)
        : FlowElement(params) 
      { 
        d_JxW = &(d_fe[0]->get_JxW());

        d_u_dphi = &(d_fe[0]->get_dphi());

        d_u_phi = &(d_fe[0]->get_phi());

        d_p_dphi = &(d_fe[1]->get_dphi());

        d_p_phi = &(d_fe[1]->get_phi());

        d_xyz = &(d_fe[0]->get_xyz());

        d_alpha_conv  = params->d_db->getDoubleWithDefault("Convection_Coefficient", 1.0);
        d_alpha_diff  = params->d_db->getDoubleWithDefault("Diffusion_Coefficient", 1.0);

      }

      /**
        Destructor.
        */
      ~NavierStokesLSWFElement() {  }

      /**
        This function is used by FlowNonlinearFEOperator to pass the address 
        of the element Input and Output vector to this class. 
        @param [in] elementInputVectors Element input vector
        @param [in] elementOutputVector Element residual vector
        */
      void setElementVectors( const std::vector<std::vector<double> > & elementInputVectors, 
          std::vector<double> & elementOutputVector )
      {
        d_elementInputVectors = elementInputVectors;
        d_elementOutputVector = &(elementOutputVector);
      }

      /**
        Element residual vector computation.
        */
      void apply() ; 

      void initTransportModel();

    protected :

      double d_density , d_fmu , d_Re ;

      const std::vector<Real> *d_JxW; 

      const std::vector<std::vector<RealGradient> > *d_u_dphi; 

      const std::vector<std::vector<Real> > *d_u_phi; 

      const std::vector<std::vector<RealGradient> > *d_p_dphi; 

      const std::vector<std::vector<Real> > *d_p_phi; 

      const std::vector<Point> *d_xyz; 

      std::vector<std::vector<double> > d_elementInputVectors; 

      std::vector<double> *d_elementOutputVector; 
   
      bool d_alpha_conv ; 
      bool d_alpha_diff ;
 
    private :

  };

}
}

#endif


