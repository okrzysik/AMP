/*
 * Steel316_MSRZC_09.h
 *
 *  Created on: Mar 11, 2010
 *      Author: bm, gad
 */

#include "AMP/materials/Steel316_MSRZC_09.h"
#include "AMP/materials/Material.h"
#include "AMP/materials/Property.h"

#include <string>
#include <vector>

namespace AMP::Materials {

namespace Steel316_MSRZC_09_NS {

//=================== Constants =====================================================

static const char *source = "Bogdan Mihaila, Marius Stan, Juan Ramirez, Alek Zubelewicz, Petrica "
                            "Cristea, Journal of Nuclear Materials 394 (2009) 182--189";

static std::initializer_list<double> TCparams = { 7.956, 1.919e-2, -3.029e-6 };
static std::initializer_list<double> DEparams = { 7989, 0.127, 1.51e-5 };
static std::initializer_list<double> TEparams = { 15.046e-6, 5.082e-9, -1.014e-12 };
static std::initializer_list<double> HCparams = { 500, 0.072, -6.37e-4, 1.73e-6 };
static std::initializer_list<double> YMparams = { 2.116e11, -5.173e7, -1.928e4 };

static const double PRatio = 0.290;

static std::initializer_list<std::string> arguments = { "temperature" };

static const double TminVal                                = 299.9;
static const double TmaxVal                                = 1E6;
static std::initializer_list<std::array<double, 2>> ranges = { { TminVal, TmaxVal } };

//=================== Classes =======================================================

class ThermalConductivityProp : public Property
{
public:
    ThermalConductivityProp()
        : Property( "Steel316_MSRZC_09_ThermalConductivity", // Name string
                    source,                                  // Reference source
                    TCparams,                                // Property parameters
                    arguments,                               // Names of arguments
                    ranges                                   // Range of variables
          )
    {
    }

    double eval( const std::vector<double> &args ) override;
};

class DensityProp : public Property
{
public:
    DensityProp()
        : Property( "Steel316_MSRZC_09_Density", // Name string
                    source,                      // Reference source
                    DEparams,                    // Property parameters
                    arguments,                   // Names of arguments
                    ranges                       // Range of variables
          )
    {
    }

    double eval( const std::vector<double> &args ) override;
};

class ThermalExpansionProp : public Property
{
public:
    ThermalExpansionProp()
        : Property( "Steel316_MSRZC_09_ThermalExpansion", // Name string
                    source,                               // Reference source
                    TEparams,                             // Property parameters
                    arguments,                            // Names of arguments
                    ranges                                // Range of variables
          )
    {
    }

    double eval( const std::vector<double> &args ) override;
};

class HeatCapacityPressureProp : public Property
{
public:
    HeatCapacityPressureProp()
        : Property( "Steel316_MSRZC_09_HeatCapacityPressure", // Name string
                    source,                                   // Reference source
                    HCparams,                                 // Property parameters
                    arguments,                                // Names of arguments
                    ranges                                    // Range of variables
          )
    {
    }

    double eval( const std::vector<double> &args ) override;
};

class YoungsModulusProp : public Property
{
public:
    YoungsModulusProp()
        : Property( "Steel316_MSRZC_09_YoungsModulus", // Name string
                    source,                            // Reference source
                    YMparams,                          // Property parameters
                    arguments,                         // Names of arguments
                    ranges                             // Range of variables
          )
    {
    }

    double eval( const std::vector<double> &args ) override;
};

class PoissonRatioProp : public Property
{
public:
    PoissonRatioProp()
        : Property( "Steel316_MSRZC_09_PoissonRatio", // Name string
                    source,                           // Reference source
                    { PRatio },                       // Property parameters
                    arguments,                        // Names of arguments
                    ranges                            // Range of variables
          )
    {
    }

    double eval( const std::vector<double> &args ) override;
};

//=================== Functions =====================================================

inline double ThermalConductivityProp::eval( const std::vector<double> &args )
{
    double T              = args[0];
    std::vector<double> p = get_parameters();
    AMP_ASSERT( T > TminVal );
    double thcond = p[0] + p[1] * T + p[2] * T * T;
    return thcond;
}

inline double DensityProp::eval( const std::vector<double> &args )
{
    double T              = args[0];
    std::vector<double> p = get_parameters();
    AMP_ASSERT( T > TminVal );
    double dens = p[0] + p[1] * T + p[2] * T * T;
    return dens;
}

inline double ThermalExpansionProp::eval( const std::vector<double> &args )
{
    double T              = args[0];
    std::vector<double> p = get_parameters();
    AMP_ASSERT( T > TminVal );
    double alpha = p[0] + p[1] * T + p[2] * T * T;
    return alpha;
}

inline double HeatCapacityPressureProp::eval( const std::vector<double> &args )
{
    double T              = args[0];
    std::vector<double> p = get_parameters();
    AMP_ASSERT( T > TminVal );
    double T2     = T * T;
    double heatcp = p[0] + p[1] * T + p[2] * T2 + p[3] * T2 * T;
    return heatcp;
}

inline double YoungsModulusProp::eval( const std::vector<double> &args )
{
    double T              = args[0];
    std::vector<double> p = get_parameters();
    AMP_ASSERT( T > TminVal );
    double youngs = p[0] + p[1] * T + p[2] * T * T;
    return youngs;
}

inline double PoissonRatioProp::eval( const std::vector<double> & ) { return PRatio; }
} // namespace Steel316_MSRZC_09_NS

//=================== Materials =====================================================
// clang-format off
Steel316_MSRZC_09::Steel316_MSRZC_09()
{
    d_propertyMap["ThermalConductivity"]  = std::make_shared<Steel316_MSRZC_09_NS::ThermalConductivityProp>();
    d_propertyMap["Density"]              = std::make_shared<Steel316_MSRZC_09_NS::DensityProp>();
    d_propertyMap["HeatCapacityPressure"] = std::make_shared<Steel316_MSRZC_09_NS::HeatCapacityPressureProp>();
    d_propertyMap["ThermalExpansion"]     = std::make_shared<Steel316_MSRZC_09_NS::ThermalExpansionProp>();
    d_propertyMap["YoungsModulus"]        = std::make_shared<Steel316_MSRZC_09_NS::YoungsModulusProp>();
    d_propertyMap["PoissonRatio"]         = std::make_shared<Steel316_MSRZC_09_NS::PoissonRatioProp>();
}
// clang-format on

} // namespace AMP::Materials
