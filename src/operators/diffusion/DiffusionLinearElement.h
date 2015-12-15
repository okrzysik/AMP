#ifndef included_AMP_DiffusionLinearElement
#define included_AMP_DiffusionLinearElement

#include <vector>

#include "utils/shared_ptr.h"

/* AMP files */
#include "operators/diffusion/DiffusionElement.h"

namespace AMP {
namespace Operator {

class DiffusionLinearElement: public DiffusionElement {
public:

    explicit DiffusionLinearElement(const AMP::shared_ptr<ElementOperationParameters>& params) :
        DiffusionElement(params),
        d_elementStiffnessMatrix(NULL)
    {
        d_num_dofs = 0;
        d_transportAtGauss = params->d_db->getBoolWithDefault(
                "TransportAtGaussPoints", true);
    }

    virtual ~DiffusionLinearElement() {
    }

    void setElementStiffnessMatrix(
            std::vector<std::vector<double> > & elementStiffnessMatrix)
    {
        d_elementStiffnessMatrix = &(elementStiffnessMatrix);
    }

    void setElementVectors(
            unsigned int num_dofs,
            const std::vector<double>& localTemp,
            const std::vector<double>& localConc,
            const std::vector<double>& localBurn)
    {
        d_num_dofs = num_dofs;
        d_LocalTemperature = localTemp;
        d_LocalConcentration = localConc;
        d_LocalBurnup = localBurn;
    }

    void apply();

protected:

    std::vector<std::vector<double> > *d_elementStiffnessMatrix;

    bool d_transportAtGauss;

    std::vector<double> d_LocalTemperature;
    std::vector<double> d_LocalConcentration;
    std::vector<double> d_LocalBurnup;

private:

    unsigned int d_num_dofs;

};

}
}

#endif

