
#ifndef included_AMP_NavierStokesGalWFElement
#define included_AMP_NavierStokesGalWFElement

#include <vector>

#include "boost/shared_ptr.hpp"

/* AMP files */
#include "FlowElement.h"
#include "NavierStokesConstants.h"

namespace AMP {
namespace Operator {

  class NavierStokesGalWFElement : public FlowElement 
  {
    public :

      /**
        Constructor.
        */
      NavierStokesGalWFElement(const boost::shared_ptr<ElementOperationParameters>& params)
        : FlowElement(params) 
      { 
        d_JxW = &(d_fe->get_JxW());

        d_dphi = &(d_fe->get_dphi());

        d_phi = &(d_fe->get_phi());

        d_xyz = &(d_fe->get_xyz());

        d_alpha_conv  = params->d_db->getDoubleWithDefault("Convection_Coefficient", 1.0);
        d_alpha_diff  = params->d_db->getDoubleWithDefault("Diffusion_Coefficient", 1.0);

      }

      /**
        Destructor.
        */
      ~NavierStokesGalWFElement() {  }

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

      double d_density;

      double d_fmu;

      const std::vector<Real> *d_JxW; 

      const std::vector<std::vector<RealGradient> > *d_dphi; 

      const std::vector<std::vector<Real> > *d_phi; 

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


