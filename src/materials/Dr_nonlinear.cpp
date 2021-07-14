/*
 * Dr_nonlinear.h
 *
 *  Created on: June 24, 2010
 *	  Author: bm
 */

#include "AMP/materials/Dr_nonlinear.h"
#include "AMP/materials/Material.h"
#include "AMP/materials/Property.h"

#include <string>
#include <vector>

namespace AMP {
namespace Materials {

namespace Dr_nonlinear_NS {

//=================== Constants =====================================================


static const double TminVal = 299.9;
static const double TmaxVal =
    1E6; // DEBUG: This value was not provided, but is needed for ranges - set arbitrarily high
static const double uminVal = 0.0;
static const double umaxVal =
    1E6; // DEBUG: This value was not provided, but is needed for ranges - set arbitrarily high

//=================== Classes =======================================================

class FickCoefficientProp : public Property
{
public:
    FickCoefficientProp()
        : Property( "Dr_nonlinear_FickCoefficient", // Name string
                    "manufactured solution -- nonlinear D(r) \n"
                    "   u(r) = 1 - r3 \n"
                    "   D(r) = D0 exp( - gamma u(r) ) \n"
                    "   S(r) = - D(r) (u''(r) + u'(r)/r - gamma [u'(r)]^2) ", // Reference source
                    { -0.1, 2 },                                              // Property parameters
                    { "temperature", "concentration" }, // Names of arguments to the eval function
                    { { TminVal, TmaxVal }, { uminVal, umaxVal } } )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;
};

//=================== Functions =====================================================

inline double FickCoefficientProp::eval( std::vector<double> &args )
{
    double T = args[0];
    double u = args[1];

    const auto &p = get_parameters();
    AMP_ASSERT( T > TminVal );
    AMP_ASSERT( u >= uminVal );

    double fick = p[0] * exp( -p[1] * u );
    return fick;
}
} // namespace Dr_nonlinear_NS

//=================== Materials =====================================================

Dr_nonlinear::Dr_nonlinear()
{
    d_propertyMap = new std::map<std::string, std::shared_ptr<Property>>();
    INSERT_PROPERTY_IN_MAP( FickCoefficient, Dr_nonlinear_NS );
}
} // namespace Materials
} // namespace AMP
