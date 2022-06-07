#include "AMP/materials/Material.h"
#include "AMP/materials/MaterialList.h"
#include "AMP/materials/Property.h"

#include <string>
#include <vector>


namespace AMP::Materials {


class Dr_nonlinear_FickCoefficientProp : public Property
{
public:
    Dr_nonlinear_FickCoefficientProp()
        : Property( "Dr_nonlinear_FickCoefficient",
                    { 1 },
                    Units(),
                    "manufactured solution -- nonlinear D(r) \n"
                    "   u(r) = 1 - r3 \n"
                    "   D(r) = D0 exp( - gamma u(r) ) \n"
                    "   S(r) = - D(r) (u''(r) + u'(r)/r - gamma [u'(r)]^2) ",
                    { "temperature", "concentration" },
                    { { 299.9, 1E6 }, { 0, 1E6 } } )
    {
    }

    void eval( AMP::Array<double> &result, const AMP::Array<double> &args ) const override
    {
        const double p[] = { -0.1, 2 };
        for ( size_t i = 0; i < result.length(); i++ ) {
            double u    = args( 1, i );
            result( i ) = p[0] * exp( -p[1] * u );
        }
    }

private:
    std::vector<double> d_params;
};


Dr_nonlinear::Dr_nonlinear()
{
    d_propertyMap["FickCoefficient"] = std::make_shared<Dr_nonlinear_FickCoefficientProp>();
}


} // namespace AMP::Materials
