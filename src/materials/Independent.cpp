#include "AMP/materials/Material.h"
#include "AMP/materials/MaterialList.h"
#include "AMP/materials/Property.h"
#include "AMP/materials/ScalarProperty.h"

#include <string>


namespace AMP::Materials {


Independent::Independent()
{
    const char *source      = "none; all ones.";
    const double thermalval = 1.0;
    const double fickval    = 1.0;
    const double soretval   = 1.0;
    const double densval    = 1.0;
    const double alphaval   = 1.0;
    const double heatcpval  = 1.0;
    const double youngsval  = 1.0;
    const double pratioval  = 0.290;
    addScalarProperty( "Density", densval, {}, source );
    addScalarProperty( "ThermalConductivity", thermalval, {}, source );
    addScalarProperty( "FickCoefficient", fickval, {}, source );
    addScalarProperty( "SoretCoefficient", soretval, {}, source );
    addScalarProperty( "DTThermalConductivity", 0, {}, source );
    addScalarProperty( "DTFickCoefficient", 0, {}, source );
    addScalarProperty( "DTSoretCoefficient", 0, {}, source );
    addScalarProperty( "DxThermalConductivity", 0, {}, source );
    addScalarProperty( "DxSoretCoefficient", 0, {}, source );
    addScalarProperty( "DxFickCoefficient", 0, {}, source );
    addScalarProperty( "HeatCapacityPressure", heatcpval, {}, source );
    addScalarProperty( "ThermalExpansion", alphaval, {}, source );
    addScalarProperty( "YoungsModulus", youngsval, {}, source );
    addScalarProperty( "PoissonRatio", pratioval, {}, source );
    addScalarProperty( "ThermalDiffusionCoefficient", fickval * soretval, {}, source );
    addScalarProperty( "VectorFickCoefficient", 1.0, {}, source );
    addScalarProperty( "TensorFickCoefficient", Array<double>( { 1, 1 }, &fickval ), {}, source );
}


} // namespace AMP::Materials
