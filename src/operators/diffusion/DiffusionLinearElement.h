#ifndef included_AMP_DiffusionLinearElement
#define included_AMP_DiffusionLinearElement

#include <vector>

#include <memory>

/* AMP files */
#include "AMP/operators/diffusion/DiffusionElement.h"

namespace AMP::Operator {

class DiffusionLinearElement : public DiffusionElement
{
public:
    explicit DiffusionLinearElement( std::shared_ptr<const ElementOperationParameters> params )
        : DiffusionElement( params ), d_elementStiffnessMatrix( nullptr )
    {
        d_transportAtGauss = params->d_db->getWithDefault<bool>( "TransportAtGaussPoints", true );
    }

    virtual ~DiffusionLinearElement() {}

    void setElementStiffnessMatrix( std::vector<std::vector<double>> &elementStiffnessMatrix )
    {
        d_elementStiffnessMatrix = &elementStiffnessMatrix;
    }

    void setElementVectors( std::map<std::string, std::vector<double>> vecs )
    {
        d_localVecs = std::move( vecs );
    }

    void apply() override;

protected:
    std::vector<std::vector<double>> *d_elementStiffnessMatrix;

    bool d_transportAtGauss;

    std::map<std::string, std::vector<double>> d_localVecs;
};


} // namespace AMP::Operator

#endif
