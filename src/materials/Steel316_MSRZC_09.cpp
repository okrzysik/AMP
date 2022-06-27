#include "AMP/materials/Material.h"
#include "AMP/materials/MaterialList.h"
#include "AMP/materials/Property.h"

#include <string>
#include <vector>


namespace AMP::Materials {


Steel316_MSRZC_09::Steel316_MSRZC_09()
{
    const char *source = "Bogdan Mihaila, Marius Stan, Juan Ramirez, Alek Zubelewicz, Petrica "
                         "Cristea, Journal of Nuclear Materials 394 (2009) 182--189";
    std::vector<double> TCparams             = { 7.956, 1.919e-2, -3.029e-6 };
    std::vector<double> DEparams             = { 7989, 0.127, 1.51e-5 };
    std::vector<double> TEparams             = { 15.046e-6, 5.082e-9, -1.014e-12 };
    std::vector<double> HCparams             = { 500, 0.072, -6.37e-4, 1.73e-6 };
    std::vector<double> YMparams             = { 2.116e11, -5.173e7, -1.928e4 };
    std::vector<double> PRparams             = { 0.290, 0 };
    std::vector<std::string> args            = { "temperature" };
    std::vector<std::array<double, 2>> range = { { 299.9, 1e6 } };
    addPolynomialProperty( "ThermalConductivity", source, {}, TCparams, args, range );
    addPolynomialProperty( "Density", source, {}, DEparams, args, range );
    addPolynomialProperty( "HeatCapacityPressure", source, {}, HCparams, args, range );
    addPolynomialProperty( "ThermalExpansion", source, {}, TEparams, args, range );
    addPolynomialProperty( "YoungsModulus", source, {}, YMparams, args, range );
    addPolynomialProperty( "PoissonRatio", source, {}, PRparams, args, range );
}


} // namespace AMP::Materials
