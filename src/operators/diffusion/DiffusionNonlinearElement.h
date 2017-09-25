#ifndef included_AMP_DiffusionNonlinearElement
#define included_AMP_DiffusionNonlinearElement

#include <vector>

#include "utils/shared_ptr.h"

/* AMP files*/

#include "operators/diffusion/DiffusionConstants.h"
#include "operators/diffusion/DiffusionElement.h"

namespace AMP {
namespace Operator {

class DiffusionNonlinearElement : public DiffusionElement
{
public:
    explicit DiffusionNonlinearElement( const AMP::shared_ptr<ElementOperationParameters> &params )
        : DiffusionElement( params ),
          d_elementOutputVector( nullptr ),
          d_transportOutputVector( nullptr ),
          d_PrincipalVariable( 0 )
    {
        d_JxW = &( d_fe->get_JxW() );

        d_dphi = &( d_fe->get_dphi() );

        d_transportAtGauss = params->d_db->getBoolWithDefault( "TransportAtGaussPoints", true );
    }

    virtual ~DiffusionNonlinearElement() {}

    void setElementInputVector( const std::vector<std::vector<double>> &elementInputVectors )
    {
        d_elementInputVectors = elementInputVectors;
    }

    void setElementVectors( const std::vector<std::vector<double>> &elementInputVectors,
                            std::vector<double> &elementOutputVector )
    {
        d_elementInputVectors = elementInputVectors;
        d_elementOutputVector = &( elementOutputVector );
    }

    void setElementTransport( const std::vector<std::vector<double>> &elementInputVectors,
                              std::vector<double> &elementOutputVector )
    {
        d_elementInputVectors   = elementInputVectors;
        d_transportOutputVector = &( elementOutputVector );
    }

    void apply() override;

    void initTransportModel();

    void setPrincipalVariable( const unsigned int var ) { d_PrincipalVariable = var; }

    bool getTransportAtGauss() { return d_transportAtGauss; }

protected:
    std::vector<std::vector<double>> d_elementInputVectors;

    std::vector<double> *d_elementOutputVector;

    std::vector<double> *d_transportOutputVector;

    std::vector<std::vector<double>> d_elementOtherVectors;

    bool d_transportAtGauss;

    unsigned int d_PrincipalVariable;

private:
};
}
}

#endif
