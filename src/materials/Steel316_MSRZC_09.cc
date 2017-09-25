/*
 * Steel316_MSRZC_09.h
 *
 *  Created on: Mar 11, 2010
 *      Author: bm, gad
 */

#include "materials/Steel316_MSRZC_09.h"
#include "materials/Material.h"
#include "materials/Property.h"

#include <string>
#include <valarray>

namespace AMP {
namespace Materials {

namespace Steel316_MSRZC_09_NS {

//=================== Constants =====================================================

static const std::string name_base( "Steel316_MSRZC_09" );
static const std::string source( "\
Bogdan Mihaila, Marius Stan, Juan Ramirez, \
Alek Zubelewicz, Petrica Cristea, \
Journal of Nuclear Materials 394 (2009) 182--189" );

static const double TCparams[3] = { 7.956, 1.919e-2, -3.029e-6 };

static const double DEparams[3] = { 7989, 0.127, 1.51e-5 };
static const double TEparams[3] = { 15.046e-6, 5.082e-9, -1.014e-12 };
static const double HCparams[4] = { 500, 0.072, -6.37e-4, 1.73e-6 };

static const double YMparams[3] = { 2.116e11, -5.173e7, -1.928e4 };
static const double PRatio      = 0.290;

static const std::string arguments[1] = { "temperature" };
static const unsigned int narguments  = 1;

static const double TminVal = 299.9;
static const double TmaxVal =
    1E6; // DEBUG: This value was not provided, but is needed for ranges - set arbitrarily high

static const double ranges[1][2] = { { TminVal, TmaxVal } };

//=================== Classes =======================================================

class ThermalConductivityProp : public Property<double>
{
public:
    ThermalConductivityProp()
        : Property<double>( name_base + "_" + "ThermalConductivity", // Name string
                            source,                                  // Reference source
                            TCparams,                                // Property parameters
                            3U,                                      // Number of parameters
                            arguments,  // Names of arguments to the eval function
                            narguments, // Number of arguments
                            ranges )
    {
    } // Range of variables

    virtual double eval( std::vector<double> &args ) override;
};

class DensityProp : public Property<double>
{
public:
    DensityProp()
        : Property<double>( name_base + "_" + "Density", // Name string
                            source,                      // Reference source
                            DEparams,                    // Property parameters
                            3U,                          // Number of parameters
                            arguments,                   // Names of arguments to the eval function
                            narguments,                  // Number of arguments
                            ranges )
    {
    } // Range of variables

    virtual double eval( std::vector<double> &args ) override;
};

class ThermalExpansionProp : public Property<double>
{
public:
    ThermalExpansionProp()
        : Property<double>( name_base + "_" + "ThermalExpansion", // Name string
                            source,                               // Reference source
                            TEparams,                             // Property parameters
                            3U,                                   // Number of parameters
                            arguments,  // Names of arguments to the eval function
                            narguments, // Number of arguments
                            ranges )
    {
    } // Range of variables

    virtual double eval( std::vector<double> &args ) override;
};

class HeatCapacityPressureProp : public Property<double>
{
public:
    HeatCapacityPressureProp()
        : Property<double>( name_base + "_" + "HeatCapacityPressure", // Name string
                            source,                                   // Reference source
                            HCparams,                                 // Property parameters
                            4U,                                       // Number of parameters
                            arguments,  // Names of arguments to the eval function
                            narguments, // Number of arguments
                            ranges )
    {
    } // Range of variables

    virtual double eval( std::vector<double> &args ) override;
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

    virtual double eval( std::vector<double> &args ) override;
};

class PoissonRatioProp : public Property<double>
{
public:
    PoissonRatioProp()
        : Property<double>( name_base + "_" + "PoissonRatio", // Name string
                            source,                           // Reference source
                            &PRatio,                          // Property parameters
                            1U,                               // Number of parameters
                            arguments, // Names of arguments to the eval function
                            0 )
    {
    } // Number of arguments

    virtual double eval( std::vector<double> &args ) override;
};

//=================== Functions =====================================================

inline double ThermalConductivityProp::eval( std::vector<double> &args )
{
    double T                = args[0];
    std::valarray<double> p = get_parameters();
    AMP_ASSERT( T > TminVal );
    double thcond = p[0] + p[1] * T + p[2] * T * T;
    return thcond;
}

inline double DensityProp::eval( std::vector<double> &args )
{
    double T                = args[0];
    std::valarray<double> p = get_parameters();
    AMP_ASSERT( T > TminVal );
    double dens = p[0] + p[1] * T + p[2] * T * T;
    return dens;
}

inline double ThermalExpansionProp::eval( std::vector<double> &args )
{
    double T                = args[0];
    std::valarray<double> p = get_parameters();
    AMP_ASSERT( T > TminVal );
    double alpha = p[0] + p[1] * T + p[2] * T * T;
    return alpha;
}

inline double HeatCapacityPressureProp::eval( std::vector<double> &args )
{
    double T                = args[0];
    std::valarray<double> p = get_parameters();
    AMP_ASSERT( T > TminVal );
    double T2     = T * T;
    double heatcp = p[0] + p[1] * T + p[2] * T2 + p[3] * T2 * T;
    return heatcp;
}

inline double YoungsModulusProp::eval( std::vector<double> &args )
{
    double T                = args[0];
    std::valarray<double> p = get_parameters();
    AMP_ASSERT( T > TminVal );
    double youngs = p[0] + p[1] * T + p[2] * T * T;
    return youngs;
}

inline double PoissonRatioProp::eval( std::vector<double> & ) { return PRatio; }
} // namespace Steel316_MSRZC_09_NS

//=================== Materials =====================================================

Steel316_MSRZC_09::Steel316_MSRZC_09()
{
    d_propertyMap = new std::map<std::string, AMP::shared_ptr<Property<double>>>();
    INSERT_PROPERTY_IN_MAP( ThermalConductivity, Steel316_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( Density, Steel316_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( HeatCapacityPressure, Steel316_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( ThermalExpansion, Steel316_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( YoungsModulus, Steel316_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( PoissonRatio, Steel316_MSRZC_09_NS );
}
} // namespace Materials
} // namespace AMP
