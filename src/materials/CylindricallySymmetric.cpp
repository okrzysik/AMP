#include "AMP/materials/CylindricallySymmetric.h"
#include "AMP/materials/Material.h"
#include "AMP/materials/MaterialList.h"
#include "AMP/materials/Property.h"
#include "AMP/materials/TensorProperty.h"
#include "AMP/materials/VectorProperty.h"
#include "AMP/utils/Utilities.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>


static inline int round_zero( double x ) { return x >= 0.0 ? ceil( x ) : floor( x ); }

namespace AMP::Materials {

/**
 * This is a nearly generic analytic model for the cylindrically symmetric diffusion in 3D.
 * It is primarily meant for manufactured solutions.
 *
 * For radial FickCoefficient, specify n+1 parameters {p0, ..., pn}. These will
 * be the coefficients of the powers of r from 0 to n.
 *
 * For longitudinal FickCoefficient, specify separate n+1 parameters {p0, ..., pn}. These will
 * be the coefficients of the powers of z from 0 to n.
 *
 * The n's for radial and longitudinal can be different.
 *
 * A tensor diffusion coefficient of the form
 *
 *	   \f$
       \left[ \begin{array}{ c c }
       k_r(r) \cos\theta^2 & k_r(r) \sin\theta \cos\theta & 0 \\
       k_r(r) \sin\theta \cos\theta & k_r(r) \sin\theta^2 & 0 \\
       0 & 0 & k_z(z)
       \end{array} \right]
       \f$
 *
 * is generated. The only thing more general is if \f$k_r\f$ and \f$k_z\f$ are functions
 * of both \f$r\f$ and \f$z\f$.
 *
 */


// Constants
static const double rMinVal = 0.0;
static const double rMaxVal = std::numeric_limits<double>::max();
static const double tMinVal = 0.0;
static const double tMaxVal = 6.283185307179586; // 2*pi
static const double zMinVal = -std::numeric_limits<double>::max();
static const double zMaxVal = std::numeric_limits<double>::max();


// Evaluate a polynomial and it's deriviative
static std::array<double, 2> evalPoly( const std::vector<double> &p, double x )
{
    AMP_ASSERT( !p.empty() );
    double y = p.back();
    for ( size_t i = p.size() - 1; i > 0; i-- )
        y = y * x + p[i - 1];
    double d = ( p.size() - 1 ) * p.back();
    for ( size_t i = p.size() - 1; i > 1; i-- )
        d = d * x + ( i - 1 ) * p[i - 1];
    return { y, d };
}


// full cylindrically symmetric tensor diffusion coefficient
CylindricallySymmetricTensor::CylindricallySymmetricTensor( const std::string &name,
                                                            std::vector<double> params )
    : TensorProperty( name,
                      "",
                      { "radius", "theta", "zee" },
                      { { rMinVal, rMaxVal }, { tMinVal, tMaxVal }, { zMinVal, zMaxVal } },
                      { 3, 3 } )
{
    AMP_ASSERT( params.size() >= 3 );
    d_variableDimensions = true;
    d_AuxiliaryDataInteger.insert( std::make_pair( "derivative", 0 ) );
    size_t nRadial = round_zero( params[0] );
    AMP_ASSERT( nRadial < params.size() - 1 );
    size_t nLongitudinal = params.size() - 1 - nRadial;
    d_paramsRadial       = std::vector<double>( &params[1], &params[1] + nRadial );
    d_paramsLongitudinal =
        std::vector<double>( &params[1 + nRadial], &params[1 + nRadial] + nLongitudinal );
}
std::vector<std::vector<double>>
CylindricallySymmetricTensor::evalTensor( const std::vector<double> &args ) const
{
    AMP_ASSERT( args.size() >= 3 );
    std::vector<std::vector<double>> result( 3, std::vector<double>( 3, 0. ) );
    auto Kr    = evalPoly( d_paramsRadial, args[0] );
    auto Kz    = evalPoly( d_paramsLongitudinal, args[2] );
    double cth = cos( args[1] );
    double sth = sin( args[1] );
    int deriv  = d_AuxiliaryDataInteger.find( "derivative" )->second;
    switch ( deriv ) {
    case 0: {
        result[0][0] = cth * cth * Kr[0];
        result[0][1] = sth * cth * Kr[0];
        result[1][0] = sth * cth * Kr[0];
        result[1][1] = sth * sth * Kr[0];
        result[2][2] = Kz[0];
        break;
    }
    case 1: {
        result[0][0] = cth * cth * Kr[1];
        result[0][1] = sth * cth * Kr[1];
        result[1][0] = sth * cth * Kr[1];
        result[1][1] = sth * sth * Kr[1];
        result[2][2] = 0.;
        break;
    }
    case 2: {
        result[0][0] = 0.;
        result[0][1] = 0.;
        result[1][0] = 0.;
        result[1][1] = 0.;
        result[2][2] = Kz[1];
        break;
    }
    default:
        AMP_ASSERT( false );
        break;
    }
    return result;
}


//=================== Materials =====================================================
CylindricallySymmetric::CylindricallySymmetric()
{
    addPolynomialProperty(
        "ScalarRadialFick", "", {}, { 1.0, 0.0 }, { "radius" }, { { rMinVal, rMaxVal } } );
    addProperty<CylindricallySymmetricTensor>( "TensorFick" );
}


} // namespace AMP::Materials
