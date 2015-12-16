//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   materials/FixedFuel.cc
 * \author Aaron Phillippe
 * \brief  Implementation file for constant fuel properties
 */
//---------------------------------------------------------------------------//

#include "FixedFuel.h"

#include "Material.h"
#include "Property.h"
#include "TensorProperty.h"
#include "VectorProperty.h"

#include <string>

namespace AMP {
namespace Materials {

namespace FixedFuel_NS {

//  =================== Constants =====================================================

static const std::string name_base( "FixedFuel" );
static const std::string source( "average values from matpro; as defined by Phillippe." );

static const double thermalval = 3.3;
static const double fickval    = 1.;
static const double soretval   = 1.;

static const double densval   = 10540.;
static const double alphaval  = 1.;
static const double heatcpval = 1.;

static const double youngsval = 1.;
static const double pratioval = 0.290;

static const double fickVectorVal[3]     = { 1., 1., 1. };
static const double fickTensorVal[3 * 3] = { 1., 1., 1., 1., 1., 1., 1., 1., 1. };

static const std::string arguments[1];

//  =================== Classes =======================================================

class ThermalConductivityProp : public Property<double>
{
public:
    ThermalConductivityProp()
        : Property<double>( name_base + "_" + "ThermalConductivity", // Name string
                            source,                                  // Reference source
                            &thermalval,                             // Property parameters
                            1U,                                      // Number of parameters
                            arguments, // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    virtual double eval(std::vector<double> &args) override;
};

class FickCoefficientProp : public Property<double>
{
public:
    FickCoefficientProp()
        : Property<double>( name_base + "_" + "FickCoefficient", // Name string
                            source,                              // Reference source
                            &fickval,                            // Property parameters
                            1U,                                  // Number of parameters
                            arguments, // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    virtual double eval(std::vector<double> &args) override;
};

class SoretCoefficientProp : public Property<double>
{
public:
    SoretCoefficientProp()
        : Property<double>( name_base + "_" + "SoretCoefficient", // Name string
                            source,                               // Reference source
                            &fickval,                             // Property parameters
                            1U,                                   // Number of parameters
                            arguments, // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    virtual double eval(std::vector<double> &args) override;
};

class DensityProp : public Property<double>
{
public:
    DensityProp()
        : Property<double>( name_base + "_" + "Density", // Name string
                            source,                      // Reference source
                            &densval,                    // Property parameters
                            1U,                          // Number of parameters
                            arguments,                   // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    virtual double eval(std::vector<double> &args) override;
};

class ThermalExpansionProp : public Property<double>
{
public:
    ThermalExpansionProp()
        : Property<double>( name_base + "_" + "ThermalExpansion", // Name string
                            source,                               // Reference source
                            &alphaval,                            // Property parameters
                            1U,                                   // Number of parameters
                            arguments, // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    virtual double eval(std::vector<double> &args) override;
};

class HeatCapacityPressureProp : public Property<double>
{
public:
    HeatCapacityPressureProp()
        : Property<double>( name_base + "_" + "HeatCapacityPressure", // Name string
                            source,                                   // Reference source
                            &heatcpval,                               // Property parameters
                            1U,                                       // Number of parameters
                            arguments, // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    virtual double eval(std::vector<double> &args) override;
};

class YoungsModulusProp : public Property<double>
{
public:
    YoungsModulusProp()
        : Property<double>( name_base + "_" + "YoungsModulus", // Name string
                            source,                            // Reference source
                            &youngsval,                        // Property parameters
                            1U,                                // Number of parameters
                            arguments, // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    virtual double eval(std::vector<double> &args) override;
};

class PoissonRatioProp : public Property<double>
{
public:
    PoissonRatioProp()
        : Property<double>( name_base + "_" + "PoissonRatio", // Name string
                            source,                           // Reference source
                            &pratioval,                       // Property parameters
                            1U,                               // Number of parameters
                            arguments, // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    virtual double eval(std::vector<double> &args) override;
};

class DTThermalConductivityProp : public Property<double>
{
public:
    DTThermalConductivityProp()
        : Property<double>( name_base + "_" + "DTThermalConductivity", // Name string
                            source,                                    // Reference source
                            &thermalval,                               // Property parameters
                            1U,                                        // Number of parameters
                            arguments, // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    virtual double eval(std::vector<double> &args) override;
};

class DTFickCoefficientProp : public Property<double>
{
public:
    DTFickCoefficientProp()
        : Property<double>( name_base + "_" + "DTFickCoefficient", // Name string
                            source,                                // Reference source
                            &fickval,                              // Property parameters
                            1U,                                    // Number of parameters
                            arguments, // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    virtual double eval(std::vector<double> &args) override;
};

class DTSoretCoefficientProp : public Property<double>
{
public:
    DTSoretCoefficientProp()
        : Property<double>( name_base + "_" + "DTSoretCoefficient", // Name string
                            source,                                 // Reference source
                            &soretval,                              // Property parameters
                            1U,                                     // Number of parameters
                            arguments, // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    virtual double eval(std::vector<double> &args) override;
};

class DxThermalConductivityProp : public Property<double>
{
public:
    DxThermalConductivityProp()
        : Property<double>( name_base + "_" + "DxThermalConductivity", // Name string
                            source,                                    // Reference source
                            &thermalval,                               // Property parameters
                            1U,                                        // Number of parameters
                            arguments, // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    virtual double eval(std::vector<double> &args) override;
};

class DxFickCoefficientProp : public Property<double>
{
public:
    DxFickCoefficientProp()
        : Property<double>( name_base + "_" + "DxFickCoefficient", // Name string
                            source,                                // Reference source
                            &fickval,                              // Property parameters
                            1U,                                    // Number of parameters
                            arguments, // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    virtual double eval(std::vector<double> &args) override;
};

class DxSoretCoefficientProp : public Property<double>
{
public:
    DxSoretCoefficientProp()
        : Property<double>( name_base + "_" + "DxSoretCoefficient", // Name string
                            source,                                 // Reference source
                            &soretval,                              // Property parameters
                            1U,                                     // Number of parameters
                            arguments, // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    virtual double eval(std::vector<double> &args) override;
};

class VectorFickCoefficientProp : public VectorProperty<double>
{
public:
    explicit VectorFickCoefficientProp( const size_t dim = 1 )
        : VectorProperty<double>( name_base + "_" + "VectorFickCoefficient", // Name string
                                  source,                                    // Reference source
                                  fickVectorVal,                             // Property parameters
                                  1U,                                        // Number of parameters
                                  arguments, // Names of arguments to the eval function
                                  0,         // Number of arguments
                                  ( double ( * )[2] )( NULL ), // ranges
                                  dim )                        // dimension
    {
        AMP_INSIST( d_nparams == dim, "dimensions and number of parameters don't match" );
        d_variableNumberParameters = true;
        d_variableDimension        = true;
    }

    // NOTE: must change dimension first before changing number of parameters
    virtual void
    set_parameters_and_number(const double *params,
                              const unsigned int nparams) override {
      AMP_INSIST(d_dimension == nparams,
                 "number of new parameters must be same as dimension");
      Property<double>::set_parameters_and_number(params, nparams);
    }

    virtual std::vector<double> evalVector(std::vector<double> &args) override;
};

class TensorFickCoefficientProp : public TensorProperty<double>
{
public:
    explicit TensorFickCoefficientProp( const std::vector<size_t> &dims = std::vector<size_t>( 2,
                                                                                               1 ) )
        : TensorProperty<double>( name_base + "_" + "TensorFickCoefficient", // Name string
                                  source,                                    // Reference source
                                  fickTensorVal,                             // Property parameters
                                  1U,                                        // Number of parameters
                                  arguments, // Names of arguments to the eval function
                                  0,         // Number of arguments
                                  ( double ( * )[2] )( NULL ), // ranges
                                  dims )                       // dimensions
    {
        AMP_INSIST( d_nparams == dims[0] * dims[1],
                    "dimensions and number of parameters don't match" );
        d_variableNumberParameters = true;
        d_variableDimensions       = true;
    }

    // NOTE: must change dimension first before changing number of parameters
    virtual void
    set_parameters_and_number(const double *params,
                              const unsigned int nparams) override {
      AMP_INSIST(d_dimensions[0] * d_dimensions[1] == nparams,
                 "number of new parameters must be product of dimensions");
      Property<double>::set_parameters_and_number(params, nparams);
    }

    virtual std::vector<std::vector<double>>
    evalTensor(std::vector<double> &args) override;
};

static const unsigned int numberThDiffParams             = 2;
static double thermalDiffusionParams[numberThDiffParams] = { 1., 1. };

static const std::string thermDiffArgs[1];
static const unsigned int numberThermDiffArgs = 0;
static const double thermDiffRanges[1][2]     = {};

#define THERMAL_DIFFUSION_DERIVATIVE
#include "ThermalDiffusionCoefficientProp.h"
#undef THERMAL_DIFFUSION_DERIVATIVE

//  =================== Functions =====================================================

inline double ThermalConductivityProp::eval( std::vector<double> & ) { return get_parameters()[0]; }

inline double FickCoefficientProp::eval( std::vector<double> & ) { return get_parameters()[0]; }

inline double SoretCoefficientProp::eval( std::vector<double> & ) { return get_parameters()[0]; }

inline double DTThermalConductivityProp::eval( std::vector<double> & ) { return 0.; }

inline double DxThermalConductivityProp::eval( std::vector<double> & ) { return 0.; }

inline double DTFickCoefficientProp::eval( std::vector<double> & ) { return 0.; }

inline double DxFickCoefficientProp::eval( std::vector<double> & ) { return 0.; }

inline double DTSoretCoefficientProp::eval( std::vector<double> & ) { return 0.; }

inline double DxSoretCoefficientProp::eval( std::vector<double> & ) { return 0.; }

inline double DensityProp::eval( std::vector<double> & ) { return get_parameters()[0]; }

inline double ThermalExpansionProp::eval( std::vector<double> & ) { return get_parameters()[0]; }

inline double HeatCapacityPressureProp::eval( std::vector<double> & )
{
    return get_parameters()[0];
}

inline double YoungsModulusProp::eval( std::vector<double> & ) { return get_parameters()[0]; }

inline double PoissonRatioProp::eval( std::vector<double> & ) { return get_parameters()[0]; }

std::vector<double> VectorFickCoefficientProp::evalVector( std::vector<double> & )
{
    std::vector<double> result( d_dimension );
    for ( size_t i = 0; i < d_dimension; i++ )
        result[i]  = d_params[i];
    return result;
}

std::vector<std::vector<double>> TensorFickCoefficientProp::evalTensor( std::vector<double> & )
{
    std::vector<std::vector<double>> result( d_dimensions[0],
                                             std::vector<double>( d_dimensions[1] ) );
    for ( size_t i = 0; i < d_dimensions[0]; i++ )
        for ( size_t j   = 0; j < d_dimensions[1]; j++ )
            result[i][j] = d_params[i * d_dimensions[1] + j];
    return result;
}
}

//  =================== Materials =====================================================

FixedFuel::FixedFuel()
{
    d_propertyMap = new std::map<std::string, AMP::shared_ptr<Property<double>>>();
    INSERT_PROPERTY_IN_MAP( ThermalConductivity, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( FickCoefficient, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( SoretCoefficient, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( DTThermalConductivity, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( DTFickCoefficient, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( DTSoretCoefficient, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( DxThermalConductivity, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( DxFickCoefficient, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( DxSoretCoefficient, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( Density, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( HeatCapacityPressure, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( ThermalExpansion, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( YoungsModulus, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( PoissonRatio, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( ThermalDiffusionCoefficient, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( VectorFickCoefficient, FixedFuel_NS );
    INSERT_PROPERTY_IN_MAP( TensorFickCoefficient, FixedFuel_NS );
}
}
}
