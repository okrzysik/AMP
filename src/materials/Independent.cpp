#include "AMP/materials/Material.h"
#include "AMP/materials/MaterialList.h"
#include "AMP/materials/Property.h"
#include "AMP/materials/ScalarProperty.h"
#include "AMP/materials/TensorProperty.h"
#include "AMP/materials/VectorProperty.h"

#include <string>

namespace AMP::Materials {

//  =================== Constants =====================================================

static const char *name_base = "Independent";
static const char *source    = "none; all ones.";

static const double thermalval = 1.;
static const double fickval    = 1.;
static const double soretval   = 1.;

static const double densval   = 1.;
static const double alphaval  = 1.;
static const double heatcpval = 1.;

static const double youngsval = 1.;
static const double pratioval = 0.290;


//  =================== Classes =======================================================

namespace Independent_NS {

class VectorFickCoefficientProp : public VectorProperty
{
public:
    explicit VectorFickCoefficientProp( const size_t dim = 1 )
        : VectorProperty( "Independent_VectorFickCoefficient", // Name string
                          source,                              // Reference source
                          { 1. },                              // Property parameters
                          {},                                  // Names of arguments
                          {},                                  // ranges
                          dim )                                // dimension
    {
        AMP_INSIST( d_params.size() == dim, "dimensions and number of parameters don't match" );
        d_variableNumberParameters = true;
        d_variableDimension        = true;
    }

    // NOTE: must change dimension first before changing number of parameters
    void set_parameters_and_number( std::vector<double> params ) override
    {
        AMP_INSIST( d_dimension == params.size(),
                    "number of new parameters must be same as dimension" );
        Property::set_parameters_and_number( params );
    }

    std::vector<double> evalVector( const std::vector<double> &args ) override;
};

class TensorFickCoefficientProp : public TensorProperty
{
public:
    explicit TensorFickCoefficientProp( const std::vector<size_t> &dims = std::vector<size_t>( 2,
                                                                                               1 ) )
        : TensorProperty( "FixedFuel_TensorFickCoefficient", // Name string
                          source,                            // Reference source
                          { 1. },                            // Property parameters
                          {},                                // Names of arguments
                          {},                                // ranges
                          dims )                             // dimensions
    {
        AMP_INSIST( d_params.size() == dims[0] * dims[1],
                    "dimensions and number of parameters don't match" );
        d_variableNumberParameters = true;
        d_variableDimensions       = true;
    }

    // NOTE: must change dimension first before changing number of parameters
    void set_parameters_and_number( std::vector<double> params ) override
    {
        AMP_INSIST( d_dimensions[0] * d_dimensions[1] == params.size(),
                    "number of new parameters must be product of dimensions" );
        Property::set_parameters_and_number( params );
    }

    std::vector<std::vector<double>> evalTensor( const std::vector<double> &args ) override;
};

static std::initializer_list<double> thermalDiffusionParams = { 1., 1. };
static std::initializer_list<std::string> thermDiffArgs     = {};
static std::vector<std::array<double, 2>> thermDiffRanges   = {};

#include "ThermalDiffusionCoefficientProp.h"


std::vector<double> VectorFickCoefficientProp::evalVector( const std::vector<double> & )
{
    std::vector<double> result( d_dimension );
    for ( size_t i = 0; i < d_dimension; i++ )
        result[i] = d_params[i];
    return result;
}

std::vector<std::vector<double>>
TensorFickCoefficientProp::evalTensor( const std::vector<double> & )
{
    std::vector<std::vector<double>> result( d_dimensions[0],
                                             std::vector<double>( d_dimensions[1] ) );
    for ( size_t i = 0; i < d_dimensions[0]; i++ )
        for ( size_t j = 0; j < d_dimensions[1]; j++ )
            result[i][j] = d_params[i * d_dimensions[1] + j];
    return result;
}
} // namespace Independent_NS

//  =================== Materials =====================================================
Independent::Independent()
{
    registerPolynomialProperty( "Density", source, std::vector<double>( { densval } ) );
    registerPolynomialProperty(
        "ThermalConductivity", source, std::vector<double>( { thermalval } ) );
    registerPolynomialProperty( "FickCoefficient", source, std::vector<double>( { fickval } ) );
    registerPolynomialProperty( "SoretCoefficient", source, std::vector<double>( { fickval } ) );
    registerPolynomialProperty( "DTThermalConductivity", source, std::vector<double>( { 0 } ) );
    registerPolynomialProperty( "DTFickCoefficient", source, std::vector<double>( { 0 } ) );
    registerPolynomialProperty( "DTSoretCoefficient", source, std::vector<double>( { 0 } ) );
    registerPolynomialProperty( "DxThermalConductivity", source, std::vector<double>( { 0 } ) );
    registerPolynomialProperty( "DxSoretCoefficient", source, std::vector<double>( { 0 } ) );
    registerPolynomialProperty( "DxFickCoefficient", source, std::vector<double>( { 0 } ) );
    registerPolynomialProperty(
        "HeatCapacityPressure", source, std::vector<double>( { heatcpval } ) );
    registerPolynomialProperty( "ThermalExpansion", source, std::vector<double>( { alphaval } ) );
    registerPolynomialProperty( "YoungsModulus", source, std::vector<double>( { youngsval } ) );
    registerPolynomialProperty( "PoissonRatio", source, std::vector<double>( { pratioval } ) );
    registerPolynomialProperty(
        "ThermalDiffusionCoefficient", source, std::vector<double>( { 0 } ) );
    d_propertyMap["ThermalDiffusionCoefficient"] =
        std::make_shared<Independent_NS::ThermalDiffusionCoefficientProp>(
            property( "FickCoefficient" ), property( "SoretCoefficient" ) );
    d_propertyMap["VectorFickCoefficient"] =
        std::make_shared<Independent_NS::VectorFickCoefficientProp>();
    d_propertyMap["TensorFickCoefficient"] =
        std::make_shared<Independent_NS::TensorFickCoefficientProp>();
}

} // namespace AMP::Materials
