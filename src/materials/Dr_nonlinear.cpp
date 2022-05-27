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
                    Units(),
                    "manufactured solution -- nonlinear D(r) \n"
                    "   u(r) = 1 - r3 \n"
                    "   D(r) = D0 exp( - gamma u(r) ) \n"
                    "   S(r) = - D(r) (u''(r) + u'(r)/r - gamma [u'(r)]^2) ",
                    { "temperature", "concentration" },
                    { { 299.9, 1E6 }, { 0, 1E6 } } ),
          d_params{ -0.1, 2 }
    {
    }

    double eval( const std::vector<double> &args ) override
    {
        double T = args[0];
        double u = args[1];

        const auto &p = d_params;
        AMP_ASSERT( T > 299.9 );
        AMP_ASSERT( u >= 0 );

        double fick = p[0] * exp( -p[1] * u );
        return fick;
    }

private:
    std::vector<double> d_params;
};


Dr_nonlinear::Dr_nonlinear()
{
    d_propertyMap["FickCoefficient"] = std::make_shared<Dr_nonlinear_FickCoefficientProp>();
}


} // namespace AMP::Materials
