/*
 * Ox_MSRZC_09.h
 *
 *  Created on: Mar 11, 2010
 *	  Author: bm, gad
 *
 * gives soret as fick * soret
 *
 */

#include "AMP/materials/Ox_MSRZC_09.h"
#include "AMP/materials/Material.h"
#include "AMP/materials/Property.h"

#include <cmath>
#include <string>

namespace AMP::Materials {

namespace Ox_MSRZC_09_NS {

//=================== Constants =====================================================

static const char *name_base = "Ox_MSRZC_09";
static const char *source    = "Bogdan Mihaila, Marius Stan, Juan Ramirez, Alek Zubelewicz, "
                            "Petrica Cristea, Journal of Nuclear Materials 394 (2009) 182--189";

static std::initializer_list<double> FCparams = { -9.386, -4.26e3, 1.2e-3, 7.5e-4 };
static std::initializer_list<double> SCparams = { -9.386,  -4.26e3,   1.2e-3, 7.5e-4,
                                                  -1380.8, -134435.5, 0.0261 };

static std::initializer_list<std::string> arguments = { "temperature", "concentration" };

static const double TminVal = 299.9;
static const double TmaxVal = 1400;
static const double uminVal = 0.0;
static const double umaxVal = 0.2;


//=================== Classes =======================================================

class FickCoefficientProp : public Property
{
public:
    FickCoefficientProp()
        : Property( "Ox_MSRZC_09_FickCoefficient", // Name string
                    source,                        // Reference source
                    FCparams,                      // Property parameters
                    arguments,                     // Names of arguments to the eval function
                    { { TminVal, TmaxVal }, { uminVal, umaxVal } } )
    {
    } // Range of variables

    double eval( const std::vector<double> &args ) override;
};

class SoretCoefficientProp : public Property
{
public:
    SoretCoefficientProp()
        : Property( "Ox_MSRZC_09_SoretCoefficient", // Name string
                    source,                         // Reference source
                    SCparams,                       // Property parameters
                    arguments,                      // Names of arguments to the eval function
                    { { TminVal, TmaxVal }, { uminVal, umaxVal } } )
    {
    } // Range of variables

    double eval( const std::vector<double> &args ) override;
};

//=================== Functions =====================================================

inline double FickCoefficientProp::eval( const std::vector<double> &args )
{
    double T = args[0];
    double u = args[1];

    std::vector<double> p = get_parameters();
    AMP_ASSERT( T > TminVal && T < TmaxVal );
    AMP_ASSERT( u >= uminVal && u <= umaxVal );

    double x = u;
    if ( x < 0.001 )
        x = 0.001;

    double expDC = p[0] + p[1] / T + p[2] * T * x + p[3] * T * log10( ( 2 + x ) / x );
    double fick  = exp( expDC * log( 10.0 ) );
    return fick;
}

inline double SoretCoefficientProp::eval( const std::vector<double> &args )
{
    double T = args[0];
    double u = args[1];

    std::vector<double> p = get_parameters();
    AMP_ASSERT( T > TminVal && T < TmaxVal );
    AMP_ASSERT( u >= uminVal && u <= umaxVal );

    double x = u;
    if ( x < 0.001 )
        x = 0.001;

    double Q_star = p[4] + p[5] * exp( -x / p[6] );
    double F_SC   = ( 2 + x ) / ( 2 * ( 1 - 3 * x ) * ( 1 - 2 * x ) );
    double R_ig   = 8.314;
    double soret  = x / F_SC * Q_star / ( R_ig * T * T );
    return soret;
}

//=================== Thermal Diffusion Interface ================================


static std::initializer_list<double> thermalDiffusionParams = { -9.386,  -4.26e3,   1.2e-3, 7.5e-4,
                                                                -9.386,  -4.26e3,   1.2e-3, 7.5e-4,
                                                                -1380.8, -134435.5, 0.0261 };

static std::initializer_list<std::string> thermDiffArgs   = { "temperature", "concentration" };
static std::vector<std::array<double, 2>> thermDiffRanges = { { TminVal, TmaxVal },
                                                              { uminVal, umaxVal } };

#include "ThermalDiffusionCoefficientProp.h"
} // namespace Ox_MSRZC_09_NS

//=================== Materials =====================================================

Ox_MSRZC_09::Ox_MSRZC_09()
{
    d_propertyMap = new std::map<std::string, std::shared_ptr<Property>>();
    INSERT_PROPERTY_IN_MAP( FickCoefficient, Ox_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( SoretCoefficient, Ox_MSRZC_09_NS );
    INSERT_PROPERTY_IN_MAP( ThermalDiffusionCoefficient, Ox_MSRZC_09_NS );
}
} // namespace AMP::Materials
