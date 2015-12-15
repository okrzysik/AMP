
#ifndef included_AMP_NavierStokesLSWFElement
#define included_AMP_NavierStokesLSWFElement

#include <vector>

#include "utils/shared_ptr.h"

/* AMP files */
#include "operators/flow/FlowElement.h"
#include "operators/flow/NavierStokesConstants.h"

namespace AMP {
namespace Operator {

class NavierStokesLSWFElement : public FlowElement 
{
public :

    //! Constructor.
    explicit NavierStokesLSWFElement(const AMP::shared_ptr<ElementOperationParameters>& params) :
        FlowElement(params),
        d_elementOutputVector(NULL)
    { 
        d_JxW = &(d_fe->get_JxW());
        d_dphi = &(d_fe->get_dphi());
        d_phi = &(d_fe->get_phi());
        d_xyz = &(d_fe->get_xyz());
        d_alpha_conv  = params->d_db->getDoubleWithDefault("Convection_Coefficient", 1.0);
        d_alpha_diff  = params->d_db->getDoubleWithDefault("Diffusion_Coefficient", 1.0);
        d_density = 0.;
        d_fmu     = 0.;
        d_Re      = 0.;
    }

    /**
      Destructor.
     */
    virtual ~NavierStokesLSWFElement() {  }

    /**
      This function is used by FlowNonlinearFEOperator to pass the address 
      of the element Input and Output vector to this class. 
      @param [in] elementInputVectors Element input vector
      @param [in] elementOutputVectors Element residual vector
     */
    void setElementVectors( const std::vector<double> & elementInputVectors, 
        std::vector<double> & elementOutputVectors )
    {
        d_elementInputVectors = elementInputVectors;
        d_elementOutputVector = &(elementOutputVectors);
    }

    /**
    Element residual vector computation.
    */
    void apply() ; 

    void initTransportModel();

protected :

    double d_density , d_fmu , d_Re ;

    const std::vector<Real> *d_JxW; 

    const std::vector<std::vector<RealGradient> > *d_dphi; 

    const std::vector<std::vector<Real> > *d_phi; 

    const std::vector<Point> *d_xyz; 

    std::vector<double>  d_elementInputVectors; 

    std::vector<double> *d_elementOutputVector; 

    bool d_alpha_conv ; 
    bool d_alpha_diff ;

private :

};


}
}

#endif


