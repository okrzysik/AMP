//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   materials/FixedFuel.cc
 * \author Aaron Phillippe
 * \brief  Implementation file for constant fuel properties
 */
//---------------------------------------------------------------------------//

#include "AMP/materials/Material.h"
#include "AMP/materials/MaterialList.h"
#include "AMP/materials/Property.h"
#include "AMP/materials/ScalarProperty.h"
#include "AMP/materials/TensorProperty.h"
#include "AMP/materials/VectorProperty.h"

#include <string>

namespace AMP::Materials {

//  =================== Constants =====================================================

static const char *name_base = "FixedFuel";
static const char *source    = "average values from matpro; as defined by Phillippe";

static const double thermalval = 3.3;
static const double fickval    = 1.;
static const double soretval   = 1.;

static const double densval   = 10540.;
static const double alphaval  = 1.;
static const double heatcpval = 1.;

static const double youngsval = 1.;
static const double pratioval = 0.290;


//  =================== Classes =======================================================

namespace FixedFuel_NS {

class VectorFickCoefficientProp : public VectorProperty
{
public:
    explicit VectorFickCoefficientProp( const size_t dim = 1 )
        : VectorProperty( "FixedFuel_VectorFickCoefficient", // Name string
                          source,                            // Reference source
                          { 1. },                            // Property parameters
                          {},                                // Names of arguments
                          {},                                // ranges
                          dim )                              // dimension
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
} // namespace FixedFuel_NS

//  =================== Materials =====================================================

// clang-format off
FixedFuel::FixedFuel()
{
    addPolynomialProperty( "Density", source, {}, { densval } );
    addPolynomialProperty( "ThermalConductivity", source, {}, { thermalval } );
    addPolynomialProperty( "FickCoefficient", source, {}, { fickval } );
    addPolynomialProperty( "SoretCoefficient", source, {}, { fickval } );
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
    d_propertyMap["ThermalDiffusionCoefficient"] = std::make_shared<FixedFuel_NS::ThermalDiffusionCoefficientProp>(property( "FickCoefficient" ), property( "SoretCoefficient" ));
    d_propertyMap["VectorFickCoefficient"] = std::make_shared<FixedFuel_NS::VectorFickCoefficientProp>();
    d_propertyMap["TensorFickCoefficient"] = std::make_shared<FixedFuel_NS::TensorFickCoefficientProp>();
}
// clang-format on


} // namespace AMP::Materials
