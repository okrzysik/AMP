/*
 * CylindricallySymmetric.h
 *
 *  Created on: June 11, 2010
 *	  Author: bm
 */

#include "AMP/materials/Material.h"
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


// Define and register the material
class CylindricallySymmetric : public AMP::Materials::Material
{
public:
    CylindricallySymmetric();
};
registerMaterial( CylindricallySymmetric, "CylindricallySymmetric" );


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
 * Example:
 * \code
 * double param[5] = {1., 2.3, 0., 4.5, 6.7};
 * RadialFickProp prop;
 * prop.set_parameters(param, 5);
 * std::vector args(1); args[0] = r;
 * double v = prop.eval(args); // v has the value 1. + 2.3*r + 4.5*r*r*r + 6.7*r*r*r*r
 * \endcode
 */
namespace CylindricallySymmetric_NS {

//=================== Constants =====================================================

static const double rMinVal = 0.0;
static const double rMaxVal = std::numeric_limits<double>::max();
static const double tMinVal = 0.0;
static const double tMaxVal = 6.283185307179586; // 2*pi
static const double zMinVal = -std::numeric_limits<double>::max();
static const double zMaxVal = std::numeric_limits<double>::max();


//=================== Classes =======================================================

/** radial diffusion coefficient */
class ScalarRadialFickProp : public Property
{
public:
    ScalarRadialFickProp()
        : Property( "CylindricallySymmetric_ScalarRadialFick", // Name string
                    "",                                        // Reference source
                    { 1.0 },                                   // Property parameters
                    { "radius" },              // Names of arguments to the eval function
                    { { rMinVal, rMaxVal } } ) // Ranges
    {
        d_variableNumberParameters = true;
    }

    /** returns property and derivative wrto r
     * \return [0]=property, [1]=property derivative wrto r
     */
    double eval( const std::vector<double> &args ) override;
};

/** radial diffusion coefficient */
class RadialFickProp : public VectorProperty
{
public:
    RadialFickProp()
        : VectorProperty( "CylindricallySymmetric_RadialFick", // Name string
                          "",                                  // Reference source
                          { 1.0 },                             // Property parameters
                          { "radius" },              // Names of arguments to the eval function
                          { { rMinVal, rMaxVal } } ) // Ranges
    {
        d_variableNumberParameters = true;
    }

    /** returns property and derivative wrto r
     * \return [0]=property, [1]=property derivative wrto r
     */
    std::vector<double> evalVector( const std::vector<double> &args ) override;
};

/** longitudinal diffusion coefficient */
class LongitudinalFickProp : public VectorProperty
{
public:
    LongitudinalFickProp()
        : VectorProperty( "CylindricallySymmetric_LongitudinalFick", // Name string
                          "",                                        // Reference source
                          { 1.0 },                                   // Property parameters
                          { "zee" },                 // Names of arguments to the eval function
                          { { zMinVal, zMaxVal } } ) // Ranges
    {
        d_variableNumberParameters = true;
    }

    std::vector<double> evalVector( const std::vector<double> &args ) override;
};

/** full cylindrically symmetric tensor diffusion coefficient
 *
 * The parameters are set as follows:
 * params[0] = number of parameters for radial
 * params[1]...params[ params[0] ] = parameters for radial
 * the rest are for the longitudinal
 * AuxiliaryInteger data "derivative" has values 0, 1, 2 for
 * zeroth, r- and z- derivatives, respectively.
 */
class TensorFickProp : public TensorProperty
{
public:
    TensorFickProp()
        : TensorProperty(
              "CylindricallySymmetric_TensorFick", // Name string
              "",                                  // Reference source
              { 1.0, 1.0, 1.0 },                   // Property parameters
              { "radius", "theta", "zee" },        // Names of arguments to the eval function
              { { rMinVal, rMaxVal }, { tMinVal, tMaxVal }, { zMinVal, zMaxVal } }, // ranges
              { 3, 3 } )                                                            // dimensions
    {
        d_variableNumberParameters = true;
        d_variableDimensions       = true;
        d_AuxiliaryDataInteger.insert( std::make_pair( "derivative", 0 ) );
        set_parameters_and_number( { 1.0, 1.0, 1.0 } );
    }

    // NOTE: must change dimension first before changing number of parameters
    void set_parameters_and_number( std::vector<double> params ) override
    {
        AMP_ASSERT( params.size() >= 3 );
        Property::set_parameters_and_number( params );
        d_nparamsRadial = round_zero( d_params[0] );
        AMP_ASSERT( d_nparamsRadial < params.size() - 1 );
        d_nparamsLongitudinal = params.size() - 1 - d_nparamsRadial;
        std::vector<double> paramsRadial( &d_params[1], &d_params[1] + d_nparamsRadial );
        std::vector<double> paramsLongitudinal( &d_params[1 + d_nparamsRadial],
                                                &d_params[1 + d_nparamsRadial] +
                                                    d_nparamsLongitudinal );
        d_radialK.set_parameters_and_number( paramsRadial );
        d_longitudinalK.set_parameters_and_number( paramsLongitudinal );
    }

    std::vector<std::vector<double>> evalTensor( const std::vector<double> &args ) override;

private:
    RadialFickProp d_radialK;
    LongitudinalFickProp d_longitudinalK;
    unsigned int d_nparamsRadial;
    unsigned int d_nparamsLongitudinal;
};

//=================== Functions =====================================================

inline double ScalarRadialFickProp::eval( const std::vector<double> &args )
{
    AMP_ASSERT( !args.empty() );
    double result;
    result = d_params.back();
    for ( size_t i = d_params.size() - 1; i > 0; i-- ) {
        result = result * args[0] + d_params[i - 1];
    }
    return result;
}

inline std::vector<double> RadialFickProp::evalVector( const std::vector<double> &args )
{
    AMP_ASSERT( !args.empty() );
    std::vector<double> result( 2 );
    result[0] = d_params.back();
    result[1] = ( d_params.size() - 1 ) * d_params.back();
    for ( size_t i = d_params.size() - 1; i > 0; i-- ) {
        result[0] = result[0] * args[0] + d_params[i - 1];
    }
    for ( size_t i = d_params.size() - 1; i > 1; i-- ) {
        result[1] = result[1] * args[0] + ( i - 1 ) * d_params[i - 1];
    }
    return result;
}

inline std::vector<double> LongitudinalFickProp::evalVector( const std::vector<double> &args )
{
    AMP_ASSERT( !args.empty() );
    std::vector<double> result( 2 );
    result[0] = d_params.back();
    result[1] = ( d_params.size() - 1 ) * d_params.back();
    for ( size_t i = d_params.size() - 1; i > 0; i-- ) {
        result[0] = result[0] * args[0] + d_params[i - 1];
    }
    for ( size_t i = d_params.size() - 1; i > 1; i-- ) {
        result[1] = result[1] * args[0] + ( i - 1 ) * d_params[i - 1];
    }
    return result;
}

std::vector<std::vector<double>> TensorFickProp::evalTensor( const std::vector<double> &args )
{
    AMP_ASSERT( args.size() > 2 );
    std::vector<std::vector<double>> result( 3, std::vector<double>( 3, 0. ) );
    std::vector<double> argr( 1, args[0] );
    std::vector<double> argz( 1, args[2] );
    std::vector<double> Kr = d_radialK.evalVector( argr );
    std::vector<double> Kz = d_longitudinalK.evalVector( argz );
    double cth             = cos( args[1] );
    double sth             = sin( args[1] );
    int deriv              = d_AuxiliaryDataInteger.find( "derivative" )->second;
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
} // namespace CylindricallySymmetric_NS

//=================== Materials =====================================================
// clang-format off
CylindricallySymmetric::CylindricallySymmetric()
{
    d_propertyMap["ScalarRadialFick"] = std::make_shared<CylindricallySymmetric_NS::ScalarRadialFickProp>();
    d_propertyMap["RadialFick"]       = std::make_shared<CylindricallySymmetric_NS::RadialFickProp>();
    d_propertyMap["LongitudinalFick"] = std::make_shared<CylindricallySymmetric_NS::LongitudinalFickProp>();
    d_propertyMap["TensorFick"]       = std::make_shared<CylindricallySymmetric_NS::TensorFickProp>();
}
// clang-format on

} // namespace AMP::Materials
