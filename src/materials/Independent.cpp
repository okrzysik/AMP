#include "AMP/materials/Material.h"
#include "AMP/materials/MaterialList.h"
#include "AMP/materials/Property.h"
#include "AMP/materials/ScalarProperty.h"
#include "AMP/materials/TensorProperty.h"
#include "AMP/materials/ThermalDiffusionCoefficientProp.h"
#include "AMP/materials/VectorProperty.h"

#include <string>


namespace AMP::Materials {


Independent::Independent()
{
    const char *source      = "none; all ones.";
    const double thermalval = 1.;
    const double fickval    = 1.;
    const double soretval   = 1.;
    const double densval    = 1.;
    const double alphaval   = 1.;
    const double heatcpval  = 1.;
    const double youngsval  = 1.;
    const double pratioval  = 0.290;
    addPolynomialProperty( "Density", source, {}, { densval } );
    addPolynomialProperty( "ThermalConductivity", source, {}, { thermalval } );
    addPolynomialProperty( "FickCoefficient", source, {}, { fickval } );
    addPolynomialProperty( "SoretCoefficient", source, {}, { soretval } );
    addPolynomialProperty( "DTThermalConductivity", source, {}, { 0 } );
    addPolynomialProperty( "DTFickCoefficient", source, {}, { 0 } );
    addPolynomialProperty( "DTSoretCoefficient", source, {}, { 0 } );
    addPolynomialProperty( "DxThermalConductivity", source, {}, { 0 } );
    addPolynomialProperty( "DxSoretCoefficient", source, {}, { 0 } );
    addPolynomialProperty( "DxFickCoefficient", source, {}, { 0 } );
    addPolynomialProperty( "HeatCapacityPressure", source, {}, { heatcpval } );
    addPolynomialProperty( "ThermalExpansion", source, {}, { alphaval } );
    addPolynomialProperty( "YoungsModulus", source, {}, { youngsval } );
    addPolynomialProperty( "PoissonRatio", source, {}, { pratioval } );
    std::vector<std::string> thermDiffArgs             = {};
    std::vector<std::array<double, 2>> thermDiffRanges = {};
    auto fick                                          = property( "FickCoefficient" );
    auto soret                                         = property( "SoretCoefficient" );
    addProperty<ThermalDiffusionCoefficientProp>(
        "ThermalDiffusionCoefficient", fick, soret, thermDiffArgs, thermDiffRanges );
    std::vector<size_t> dims = { 1, 1 };
    addProperty<ScalarVectorProperty>( "VectorFickCoefficient", 1.0, source );
    addProperty<ScalarTensorProperty>(
        "TensorFickCoefficient", source, dims, std::vector<double>( 1, 1.0 ) );
}


} // namespace AMP::Materials
