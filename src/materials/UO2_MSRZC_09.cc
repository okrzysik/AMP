/*
 * UO2_MSRZC_09.h
 *
 *  Created on: Mar 11, 2010
 *	  Author: bm, gad
 *
 *
 */

#include "materials/UO2_MSRZC_09.h"
#include "materials/Material.h"
#include "materials/Property.h"

#include <string>
#include <vector>

namespace AMP {
namespace Materials {

namespace UO2_MSRZC_09_NS {

//=================== Constants =====================================================

static const std::string name_base( "UO2_MSRZC_09" );
static const std::string source( "\
Bogdan Mihaila, Marius Stan, Juan Ramirez, \
Alek Zubelewicz, Petrica Cristea, \
Journal of Nuclear Materials 394 (2009) 182--189" );

static const double TCparams[5] = { 3.24e-2, 2.51e-4, 3.67, -4.73e-4, 5.95e-11 };

static const double DEparams[9] = { 0.99672,  1.179e-5,   -2.429e-9, 1.219e-12, 0.99734,
                                    9.082e-6, -2.705e-10, 4.391e-13, 10960 };
static const double TEparams[8] = { 1.1833e-5, -5.013e-9,  3.756e-12, -6.125e-17,
                                    9.828e-6,  -6.390e-10, 1.33e-12,  -1.757e-17 };
static const double HCparams[9] = { 52.174,    45.806,      87.951e-3, -7.3461e-2, -84.241e-6,
                                    31.542e-9, -2.6334e-12, -713910,   -295090 };

static const double YMparams[3] = { 2.334e11, 1.095e-4, -1.34 };
static const double PRatio      = 0.316;

static const std::string noarguments[1];
static const std::string arguments[2] = { "temperature", "concentration" };
static const unsigned int narguments  = 2;

static const double TminVal = 299.9;
static const double TmaxVal = 1400;
static const double uminVal = 0.0;
static const double umaxVal = 0.2;

static const double ranges[2][2] = { { TminVal, TmaxVal }, { uminVal, umaxVal } };

//=================== Classes =======================================================

class ThermalConductivityProp : public Property<double>
{
public:
    ThermalConductivityProp()
        : Property<double>( name_base + "_" + "ThermalConductivity", // Name string
                            source,                                  // Reference source
                            TCparams,                                // Property parameters
                            5U,                                      // Number of parameters
                            arguments,  // Names of arguments to the eval function
                            narguments, // Number of arguments
                            ranges )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;
};

class DensityProp : public Property<double>
{
public:
    DensityProp()
        : Property<double>( name_base + "_" + "Density", // Name string
                            source,                      // Reference source
                            DEparams,                    // Property parameters
                            9U,                          // Number of parameters
                            arguments,                   // Names of arguments to the eval function
                            narguments,                  // Number of arguments
                            ranges )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;
};

class ThermalExpansionProp : public Property<double>
{
public:
    ThermalExpansionProp()
        : Property<double>( name_base + "_" + "ThermalExpansion", // Name string
                            source,                               // Reference source
                            TEparams,                             // Property parameters
                            8U,                                   // Number of parameters
                            arguments,  // Names of arguments to the eval function
                            narguments, // Number of arguments
                            ranges )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;
};

class HeatCapacityPressureProp : public Property<double>
{
public:
    HeatCapacityPressureProp()
        : Property<double>( name_base + "_" + "HeatCapacityPressure", // Name string
                            source,                                   // Reference source
                            HCparams,                                 // Property parameters
                            9U,                                       // Number of parameters
                            arguments,  // Names of arguments to the eval function
                            narguments, // Number of arguments
                            ranges )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;
};

class YoungsModulusProp : public Property<double>
{
public:
    YoungsModulusProp()
        : Property<double>( name_base + "_" + "YoungsModulus", // Name string
                            source,                            // Reference source
                            YMparams,                          // Property parameters
                            3U,                                // Number of parameters
                            arguments,  // Names of arguments to the eval function
                            narguments, // Number of arguments
                            ranges )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;
};

class PoissonRatioProp : public Property<double>
{
public:
    PoissonRatioProp()
        : Property<double>( name_base + "_" + "PoissonRatio", // Name string
                            source,                           // Reference source
                            &PRatio,                          // Property parameters
                            1U,                               // Number of parameters
                            noarguments, // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    double eval( std::vector<double> &args ) override;
};

class DxThermalConductivityProp : public Property<double>
{
public:
    DxThermalConductivityProp()
        : Property<double>( name_base + "_" + "DxThermalConductivity", // Name string
                            source,                                    // Reference source
                            TCparams,                                  // Property parameters
                            5U,                                        // Number of parameters
                            arguments,  // Names of arguments to the eval function
                            narguments, // Number of arguments
                            ranges )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;
};

class DTThermalConductivityProp : public Property<double>
{
public:
    DTThermalConductivityProp()
        : Property<double>( name_base + "_" + "DTThermalConductivity", // Name string
                            source,                                    // Reference source
                            TCparams,                                  // Property parameters
                            5U,                                        // Number of parameters
                            arguments,  // Names of arguments to the eval function
                            narguments, // Number of arguments
                            ranges )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;
};


//=================== Functions =====================================================

inline double ThermalConductivityProp::eval( std::vector<double> &args )
{
    double T = args[0];
    double u = args[1];

    std::valarray<double> p = get_parameters();
    AMP_ASSERT( T > TminVal && T < TmaxVal );
    AMP_ASSERT( u >= uminVal && u <= umaxVal );

    double x = u;
    if ( x < 0.001 )
        x = 0.001;

    double lambda0 = 1. / ( p[0] + p[1] * T );
    double Dtherm  = p[2] * exp( p[3] * T );
    double theta   = Dtherm * sqrt( 2 * x * lambda0 );
    double thcond  = p[4] * T * T * T + lambda0 * atan( theta ) / theta;
    return thcond;
}

inline double DensityProp::eval( std::vector<double> &args )
{
    double T = args[0];
    double u = args[1];

    std::valarray<double> p = get_parameters();
    AMP_ASSERT( T > TminVal && T < TmaxVal );
    AMP_ASSERT( u >= uminVal && u <= umaxVal );

    double T2 = T * T;
    double expDe;
    if ( T > 923 ) {
        expDe = p[0] + p[1] * T + p[2] * T2 + p[3] * T2 * T;
    } else {
        expDe = p[4] + p[5] * T + p[6] * T2 + p[7] * T2 * T;
    }
    double dens = p[8] * exp( -3 * log( expDe ) );
    return dens;
}

inline double ThermalExpansionProp::eval( std::vector<double> &args )
{
    double T = args[0];
    double u = args[1];

    std::valarray<double> p = get_parameters();
    AMP_ASSERT( T > TminVal && T < TmaxVal );
    AMP_ASSERT( u >= uminVal && u <= umaxVal );

    double T2 = T * T;
    double alpha;
    if ( T > 923 ) {
        alpha = p[0] + p[1] * T + p[2] * T2 + p[3] * T2 * T;
    } else {
        alpha = p[4] + p[5] * T + p[6] * T2 + p[7] * T2 * T;
    }
    return alpha;
}

inline double HeatCapacityPressureProp::eval( std::vector<double> &args )
{
    double T = args[0];
    double u = args[1];

    std::valarray<double> p = get_parameters();
    AMP_ASSERT( T > TminVal && T < TmaxVal );
    AMP_ASSERT( u >= uminVal && u <= umaxVal );

    double x = u;
    if ( x < 0.001 )
        x = 0.001;

    double T2     = T * T;
    double heatcp = p[0] + p[1] * x + ( p[2] + p[3] * x ) * T +
                    ( 1 - x ) * ( p[4] * T2 + p[5] * T2 * T + p[6] * T2 * T2 ) +
                    ( p[7] + p[8] * x ) / T2;
    return heatcp;
}

inline double YoungsModulusProp::eval( std::vector<double> &args )
{
    double T = args[0];
    double u = args[1];

    std::valarray<double> p = get_parameters();
    AMP_ASSERT( T > TminVal && T < TmaxVal );
    AMP_ASSERT( u >= uminVal && u <= umaxVal );

    double x = u;
    if ( x < 0.001 )
        x = 0.001;

    double youngs = p[0] * ( 1. + p[1] * T ) * exp( p[2] * x );
    return youngs;
}

inline double PoissonRatioProp::eval( std::vector<double> & ) { return PRatio; }

inline double DxThermalConductivityProp::eval( std::vector<double> &args )
{
    double T = args[0];
    double u = args[1];

    std::valarray<double> p = get_parameters();
    AMP_ASSERT( T > TminVal && T < TmaxVal );
    AMP_ASSERT( u >= uminVal && u <= umaxVal );

    double x = u;
    if ( x < 0.001 )
        x = 0.001;

    double DxThCond =
        ( ( -0.35355339059327373 *
            atan( 1.4142135623730951 * exp( T * p[3] ) * sqrt( x / ( p[0] + T * p[1] ) ) * p[2] ) *
            sqrt( x / ( p[0] + T * p[1] ) ) ) /
              ( exp( T * p[3] ) * p[2] ) +
          ( 0.5 * x ) / ( p[0] + T * p[1] + 2.0 * exp( 2 * T * p[3] ) * x * p[2] * p[2] ) ) /
        x * x;

    return DxThCond;
}

inline double DTThermalConductivityProp::eval( std::vector<double> &args )
{
    double T = args[0];
    double u = args[1];

    std::valarray<double> p = get_parameters();
    AMP_ASSERT( T > TminVal && T < TmaxVal );
    AMP_ASSERT( u >= uminVal && u <= umaxVal );

    double x = u;
    if ( x < 0.001 )
        x = 0.001;

    double DTThCond =
        ( -0.5 * p[1] ) / ( ( p[0] + T * p[1] ) *
                            ( p[0] + T * p[1] + 2.0 * exp( 2 * T * p[3] ) * x * p[2] * p[2] ) ) +
        ( 1. * p[3] ) / ( p[0] + T * p[1] + 2.0 * exp( 2 * T * p[3] ) * x * p[2] * p[2] ) +
        ( atan( 1.4142135623730951 * exp( T * p[3] ) * sqrt( x / ( p[0] + T * p[1] ) ) * p[2] ) *
          ( -0.7071067811865475 * p[0] * p[3] +
            p[1] * ( -0.35355339059327373 - 0.7071067811865475 * T * p[3] ) ) ) /
            ( exp( T * p[3] ) * sqrt( x / ( p[0] + T * p[1] ) ) * ( p[0] + T * p[1] ) *
              ( p[0] + T * p[1] ) * p[2] ) +
        3 * T * T * p[4];

    return DTThCond;
}
} // namespace UO2_MSRZC_09_NS

//=================== Materials =====================================================

UO2_MSRZC_09::UO2_MSRZC_09()
{
    d_propertyMap = new std::map<std::string, AMP::shared_ptr<Property<double>>>();
    INSERT_PROPERTY_IN_MAP( ThermalConductivity, UO2_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( Density, UO2_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( HeatCapacityPressure, UO2_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( ThermalExpansion, UO2_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( YoungsModulus, UO2_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( PoissonRatio, UO2_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( DxThermalConductivity, UO2_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( DTThermalConductivity, UO2_MSRZC_09_NS );
}
} // namespace Materials
} // namespace AMP
