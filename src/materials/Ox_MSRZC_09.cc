/*
 * Ox_MSRZC_09.h
 *
 *  Created on: Mar 11, 2010
 *	  Author: bm, gad
 *
 * gives soret as fick * soret
 *
 */

#include "materials/Ox_MSRZC_09.h"
#include "materials/Material.h"
#include "materials/Property.h"

#include <math.h>
#include <string>
#include <valarray>

namespace AMP {
namespace Materials {

namespace Ox_MSRZC_09_NS {

//=================== Constants =====================================================

static const std::string name_base( "Ox_MSRZC_09" );
static const std::string source( "\
Bogdan Mihaila, Marius Stan, Juan Ramirez, \
Alek Zubelewicz, Petrica Cristea, \
Journal of Nuclear Materials 394 (2009) 182--189" );

static const double FCparams[4] = { -9.386, -4.26e3, 1.2e-3, 7.5e-4 };
static const double SCparams[7] = { -9.386, -4.26e3, 1.2e-3, 7.5e-4, -1380.8, -134435.5, 0.0261 };

static const std::string arguments[2] = { "temperature", "concentration" };
static const unsigned int narguments  = 2;

static const double TminVal = 299.9;
static const double TmaxVal = 1400;
static const double uminVal = 0.0;
static const double umaxVal = 0.2;

static const double ranges[2][2] = { { TminVal, TmaxVal }, { uminVal, umaxVal } };

//=================== Classes =======================================================

class FickCoefficientProp : public Property<double> {
public:
    FickCoefficientProp()
        : Property<double>( name_base + "_" + "FickCoefficient", // Name string
                            source,                              // Reference source
                            FCparams,                            // Property parameters
                            4U,                                  // Number of parameters
                            arguments,  // Names of arguments to the eval function
                            narguments, // Number of arguments
                            ranges )
    {
    } // Range of variables

    virtual double eval( std::vector<double> &args );
};

class SoretCoefficientProp : public Property<double> {
public:
    SoretCoefficientProp()
        : Property<double>( name_base + "_" + "SoretCoefficient", // Name string
                            source,                               // Reference source
                            SCparams,                             // Property parameters
                            7U,                                   // Number of parameters
                            arguments,  // Names of arguments to the eval function
                            narguments, // Number of arguments
                            ranges )
    {
    } // Range of variables

    virtual double eval( std::vector<double> &args );
};

//=================== Functions =====================================================

inline double FickCoefficientProp::eval( std::vector<double> &args )
{
    double T = args[0];
    double u = args[1];

    std::valarray<double> p = get_parameters();
    AMP_ASSERT( T > TminVal && T < TmaxVal );
    AMP_ASSERT( u >= uminVal && u <= umaxVal );

    double x           = u;
    if ( x < 0.001 ) x = 0.001;

    double expDC = p[0] + p[1] / T + p[2] * T * x + p[3] * T * log10( ( 2 + x ) / x );
    double fick  = exp( expDC * log( 10.0 ) );
    return fick;
}

inline double SoretCoefficientProp::eval( std::vector<double> &args )
{
    double T = args[0];
    double u = args[1];

    std::valarray<double> p = get_parameters();
    AMP_ASSERT( T > TminVal && T < TmaxVal );
    AMP_ASSERT( u >= uminVal && u <= umaxVal );

    double x           = u;
    if ( x < 0.001 ) x = 0.001;

    double Q_star = p[4] + p[5] * exp( -x / p[6] );
    double F_SC   = ( 2 + x ) / ( 2 * ( 1 - 3 * x ) * ( 1 - 2 * x ) );
    double R_ig   = 8.314;
    double soret  = x / F_SC * Q_star / ( R_ig * T * T );
    return soret;
}

//=================== Thermal Diffusion Interface ================================

static const unsigned int numberThDiffParams             = 11;
static double thermalDiffusionParams[numberThDiffParams] = { -9.386,  -4.26e3,   1.2e-3, 7.5e-4,
                                                             -9.386,  -4.26e3,   1.2e-3, 7.5e-4,
                                                             -1380.8, -134435.5, 0.0261 };

static const std::string thermDiffArgs[2]     = { "temperature", "concentration" };
static const unsigned int numberThermDiffArgs = 2;
static const double thermDiffRanges[2][2]     = { { TminVal, TmaxVal }, { uminVal, umaxVal } };

#include "ThermalDiffusionCoefficientProp.h"
}

//=================== Materials =====================================================

Ox_MSRZC_09::Ox_MSRZC_09()
{
    d_propertyMap = new std::map<std::string, AMP::shared_ptr<Property<double>>>();
    INSERT_PROPERTY_IN_MAP( FickCoefficient, Ox_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( SoretCoefficient, Ox_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( ThermalDiffusionCoefficient, Ox_MSRZC_09_NS );
}
}
}
