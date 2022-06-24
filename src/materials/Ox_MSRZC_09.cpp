#include "AMP/materials/Material.h"
#include "AMP/materials/MaterialList.h"
#include "AMP/materials/Property.h"
#include "AMP/materials/ScalarProperty.h"

#include <cmath>
#include <string>


namespace AMP::Materials {


//=================== Classes =======================================================

namespace Ox_MSRZC_09_NS {


static const char *source = "Bogdan Mihaila, Marius Stan, Juan Ramirez, Alek Zubelewicz, "
                            "Petrica Cristea, Journal of Nuclear Materials 394 (2009) 182--189";
static double FCparams[]                            = { -9.386, -4.26e3, 1.2e-3, 7.5e-4 };
static double SCparams[]                            = { -1380.8, -134435.5, 0.0261 };
static std::initializer_list<std::string> arguments = { "temperature", "concentration" };
static const double TminVal                         = 299.9;
static const double TmaxVal                         = 1400;
static const double uminVal                         = 0.0;
static const double umaxVal                         = 0.2;


double evalFick( double T, double u )
{
    AMP_ASSERT( T > TminVal && T < TmaxVal );
    AMP_ASSERT( u >= uminVal && u <= umaxVal );
    double x = u;
    if ( x < 0.001 )
        x = 0.001;
    double expDC = FCparams[0] + FCparams[1] / T + FCparams[2] * T * x +
                   FCparams[3] * T * log10( ( 2 + x ) / x );
    double fick = exp( expDC * log( 10.0 ) );
    return fick;
}
double evalSoret( double T, double u )
{
    AMP_ASSERT( T > TminVal && T < TmaxVal );
    AMP_ASSERT( u >= uminVal && u <= umaxVal );
    double x = u;
    if ( x < 0.001 )
        x = 0.001;
    double Q_star = SCparams[0] + SCparams[1] * exp( -x / SCparams[2] );
    double F_SC   = ( 2 + x ) / ( 2 * ( 1 - 3 * x ) * ( 1 - 2 * x ) );
    double R_ig   = 8.314;
    double soret  = x / F_SC * Q_star / ( R_ig * T * T );
    return soret;
}


class FickCoefficientProp : public Property
{
public:
    FickCoefficientProp( const std::string &name )
        : Property( name,
                    { 1 },
                    Units(),
                    source,
                    arguments,
                    { { TminVal, TmaxVal }, { uminVal, umaxVal } } )
    {
    }
    void eval( AMP::Array<double> &result, const AMP::Array<double> &args ) const override
    {
        for ( size_t i = 0; i < result.length(); i++ ) {
            double T    = args( 0, i );
            double u    = args( 1, i );
            result( i ) = evalFick( T, u );
        }
    }
};

class SoretCoefficientProp : public Property
{
public:
    SoretCoefficientProp( const std::string &name )
        : Property( name,
                    { 1 },
                    Units(),
                    source,
                    arguments,
                    { { TminVal, TmaxVal }, { uminVal, umaxVal } } )
    {
    }
    void eval( AMP::Array<double> &result, const AMP::Array<double> &args ) const override
    {
        for ( size_t i = 0; i < result.length(); i++ ) {
            double T    = args( 0, i );
            double u    = args( 1, i );
            result( i ) = evalSoret( T, u );
        }
    }
};


class ThermalDiffusionCoefficientProp : public Property
{
public:
    ThermalDiffusionCoefficientProp( const std::string &name )
        : Property( name,
                    { 1 },
                    Units(),
                    source,
                    arguments,
                    { { TminVal, TmaxVal }, { uminVal, umaxVal } } )
    {
    }
    void eval( AMP::Array<double> &result, const AMP::Array<double> &args ) const override
    {
        for ( size_t i = 0; i < result.length(); i++ ) {
            double T    = args( 0, i );
            double u    = args( 1, i );
            result( i ) = evalFick( T, u ) * evalSoret( T, u );
        }
    }
};


} // namespace Ox_MSRZC_09_NS


//=================== Materials =====================================================
Ox_MSRZC_09::Ox_MSRZC_09()
{
    addProperty<Ox_MSRZC_09_NS::FickCoefficientProp>( "FickCoefficient" );
    addProperty<Ox_MSRZC_09_NS::SoretCoefficientProp>( "SoretCoefficient" );
    addProperty<Ox_MSRZC_09_NS::ThermalDiffusionCoefficientProp>( "ThermalDiffusionCoefficient" );
}


} // namespace AMP::Materials
