/*
 * Dr_nonlinear.h
 *
 *  Created on: June 24, 2010
 *	  Author: bm
 */

#include "materials/Dr_nonlinear.h"
#include "materials/Material.h"
#include "materials/Property.h"

#include <string>
#include <valarray>

namespace AMP {
namespace Materials {

namespace Dr_nonlinear_NS {

//=================== Constants =====================================================

static const std::string name_base( "Dr_nonlinear" );
static const std::string source( "\
manufactured solution -- nonlinear D(r) \
u(r) = 1 - r3 \
D(r) = D0 exp( - gamma u(r) ) \
S(r) = - D(r) (u''(r) + u'(r)/r - gamma [u'(r)]^2) " );

static const double FCparams[2] = { -0.1, 2 };

static const std::string arguments[2] = { "temperature", "concentration" };
static const unsigned int narguments  = 2;

static const double TminVal = 299.9;
static const double TmaxVal =
    1E6; // DEBUG: This value was not provided, but is needed for ranges - set arbitrarily high
static const double uminVal = 0.0;
static const double umaxVal =
    1E6; // DEBUG: This value was not provided, but is needed for ranges - set arbitrarily high

static const double ranges[2][2] = { { TminVal, TmaxVal }, { uminVal, umaxVal } };

//=================== Classes =======================================================

class FickCoefficientProp : public Property<double>
{
public:
    FickCoefficientProp()
        : Property<double>( name_base + "_" + "FickCoefficient", // Name string
                            source,                              // Reference source
                            FCparams,                            // Property parameters
                            2U,                                  // Number of parameters
                            arguments,  // Names of arguments to the eval function
                            narguments, // Number of arguments
                            ranges )
    {
    } // Range of variables

    virtual double eval( std::vector<double> &args ) override;
};

//=================== Functions =====================================================

inline double FickCoefficientProp::eval( std::vector<double> &args )
{
    double T = args[0];
    double u = args[1];

    std::valarray<double> p = get_parameters();
    AMP_ASSERT( T > TminVal );
    AMP_ASSERT( u >= uminVal );

    double fick = p[0] * exp( -p[1] * u );
    return fick;
}
} // namespace Dr_nonlinear_NS

//=================== Materials =====================================================

Dr_nonlinear::Dr_nonlinear()
{
    d_propertyMap = new std::map<std::string, AMP::shared_ptr<Property<double>>>();
    INSERT_PROPERTY_IN_MAP( FickCoefficient, Dr_nonlinear_NS );
}
} // namespace Materials
} // namespace AMP
