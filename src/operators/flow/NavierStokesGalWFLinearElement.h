
#ifndef included_AMP_NavierStokesGalWFLinearElement
#define included_AMP_NavierStokesGalWFLinearElement

#include <vector>

#include <memory>

/* AMP files */
#include "AMP/operators/flow/FlowElement.h"
#include "AMP/operators/flow/NavierStokesConstants.h"

namespace AMP::Operator {

class NavierStokesGalWFLinearElement : public FlowElement
{
public:
    explicit NavierStokesGalWFLinearElement(
        std::shared_ptr<const ElementOperationParameters> params )
        : FlowElement( params ), d_elementStiffnessMatrix( nullptr )
    {
        d_JxW     = &( d_fe->get_JxW() );
        d_dphi    = &( d_fe->get_dphi() );
        d_phi     = &( d_fe->get_phi() );
        d_xyz     = &( d_fe->get_xyz() );
        d_density = 0.;
        d_fmu     = 0.;
    }

    virtual ~NavierStokesGalWFLinearElement() {}

    void setElementStiffnessMatrix( std::vector<std::vector<double>> &elementStiffnessMatrix )
    {
        d_elementStiffnessMatrix = &( elementStiffnessMatrix );
    }

    void setElementVectors( const std::vector<std::vector<double>> &elementInputVectors )
    {
        d_elementInputVectors = elementInputVectors;
    }

    void apply() override;

protected:
    double d_density;

    double d_fmu;

    const std::vector<libMesh::Real> *d_JxW;

    const std::vector<std::vector<libMesh::RealGradient>> *d_dphi;

    const std::vector<std::vector<libMesh::Real>> *d_phi;

    const std::vector<libMesh::Point> *d_xyz;

    std::vector<std::vector<double>> d_elementInputVectors;

    std::vector<std::vector<double>> *d_elementStiffnessMatrix;

private:
};
} // namespace AMP::Operator

#endif
