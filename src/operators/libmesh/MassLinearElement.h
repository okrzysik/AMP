#ifndef included_AMP_MassLinearElement
#define included_AMP_MassLinearElement

#include <vector>

#include "AMP/utils/shared_ptr.h"

// AMP files
#include "AMP/operators/libmesh/MassElement.h"

namespace AMP {
namespace Operator {


class MassLinearElement : public MassElement
{
public:
    explicit MassLinearElement( const AMP::shared_ptr<ElementOperationParameters> &params )
        : MassElement( params ),
          d_elementMassMatrix( NULL ),
          d_equation( MassDensityModel::MassEquation::UnknownMassEquation )
    {
        d_densityAtGauss = params->d_db->getWithDefault( "DensityAtGaussPoints", true );
    }

    virtual ~MassLinearElement() {}

    void setElementMassMatrix( std::vector<std::vector<double>> &elementMassMatrix )
    {
        d_elementMassMatrix = &( elementMassMatrix );
    }

    void setElementVectors( const std::vector<double> &localTemp,
                            const std::vector<double> &localConc,
                            const std::vector<double> &localBurn )
    {
        d_LocalTemperature   = localTemp;
        d_LocalConcentration = localConc;
        d_LocalBurnup        = localBurn;
    }

    void apply();

protected:
    std::vector<std::vector<double>> *d_elementMassMatrix;

    bool d_densityAtGauss;

    std::vector<double> d_LocalTemperature;
    std::vector<double> d_LocalConcentration;
    std::vector<double> d_LocalBurnup;

    MassDensityModel::MassEquation d_equation;

private:
};
} // namespace Operator
} // namespace AMP

#endif
