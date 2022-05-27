#include "AMP/materials/Material.h"
#include "AMP/materials/MaterialList.h"
#include "AMP/materials/Property.h"
#include "AMP/materials/ScalarProperty.h"
#include "AMP/materials/TensorProperty.h"
#include "AMP/materials/ThermalDiffusionCoefficientProp.h"
#include "AMP/materials/VectorProperty.h"
#include <cmath>
#include <string>

namespace AMP::Materials {

//=================== Constants =====================================================

static const char *source = "Bogdan Mihaila, Marius Stan, Juan Ramirez, Alek Zubelewicz, "
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

namespace Ox_MSRZC_09_NS {

class FickCoefficientProp : public Property
{
public:
    FickCoefficientProp( const std::string &name )
        : Property( name,
                    Units(),
                    source,
                    FCparams,
                    arguments,
                    { { TminVal, TmaxVal }, { uminVal, umaxVal } } )
    {
    } // Range of variables

    double eval( const std::vector<double> &args ) override
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
};

class SoretCoefficientProp : public Property
{
public:
    SoretCoefficientProp( const std::string &name )
        : Property( name,
                    Units(),
                    source,
                    SCparams,
                    arguments,
                    { { TminVal, TmaxVal }, { uminVal, umaxVal } } )
    {
    } // Range of variables

    double eval( const std::vector<double> &args ) override
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
};


} // namespace Ox_MSRZC_09_NS


//=================== Materials =====================================================
Ox_MSRZC_09::Ox_MSRZC_09()
{
    addProperty<Ox_MSRZC_09_NS::FickCoefficientProp>( "FickCoefficient" );
    addProperty<Ox_MSRZC_09_NS::SoretCoefficientProp>( "SoretCoefficient" );
    std::vector<double> thermalDiffusionParams         = { -9.386,  -4.26e3,   1.2e-3, 7.5e-4,
                                                   -9.386,  -4.26e3,   1.2e-3, 7.5e-4,
                                                   -1380.8, -134435.5, 0.0261 };
    std::vector<std::string> thermDiffArgs             = { "temperature", "concentration" };
    std::vector<std::array<double, 2>> thermDiffRanges = { { TminVal, TmaxVal },
                                                           { uminVal, umaxVal } };
    auto fick                                          = property( "FickCoefficient" );
    auto soret                                         = property( "SoretCoefficient" );
    addProperty<ThermalDiffusionCoefficientProp>( "ThermalDiffusionCoefficient",
                                                  fick,
                                                  soret,
                                                  thermalDiffusionParams,
                                                  thermDiffArgs,
                                                  thermDiffRanges );
}


} // namespace AMP::Materials
