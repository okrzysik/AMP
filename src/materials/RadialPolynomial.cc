/*
 * RadialPolynomial.h
 *
 *  Created on: June 11, 2010
 *	  Author: bm
 */

#include "RadialPolynomial.h"

#include "materials/Property.h"
#include "materials/Material.h"
#include "utils/Utilities.h"

#include <string>
#include <valarray>
#include <cmath>
#include <limits>

namespace AMP { 
namespace Materials {

/**
 * This is a generic analytic model for the cylindrically symmetric diffusion in 3D.
 * For each property in this file, specify n+1 parameters {p0, ..., pn}. These will
 * be the coefficients of the powers of r from 0 to n.
 *
 * Example:
 * \code
 * double param[5] = {1., 2.3, 0., 4.5, 6.7};
 * FickCoefficientProp prop;
 * prop.set_parameters(param, 5);
 * std::vector args(1); args[0] = r;
 * double v = prop.eval(args); // v has the value 1. + 2.3*r + 4.5*r*r*r + 6.7*r*r*r*r
 * \endcode
 */
namespace RadialPolynomial_NS {

//=================== Constants =====================================================

	static const std::string name_base("RadialPolynomial");
	static const std::string source("");

	static const double rfMinVal=0.0;
	static const double rfMaxVal=std::numeric_limits<double>::max();

	static const std::string argumentsFick[1]={"radius"};
	static const unsigned int nargumentsFick = 1;
	static const double fickRanges[1][2]={{rfMinVal, rfMaxVal}};

	static const size_t nparams=2;
	static const double params[] = {0., 1.};

//=================== Classes =======================================================

	class FickCoefficientProp : public AMP::Materials::Property<double>{
	public:
		FickCoefficientProp() :
		  AMP::Materials::Property<double> (	name_base + "_" + "FickCoefficient",		// Name string
								source,										// Reference source
								params,										// Property parameters
								nparams,									// Number of parameters
								argumentsFick,								// Names of arguments to the eval function
								nargumentsFick,								// Number of arguments
								fickRanges)									// Ranges
								{d_variableNumberParameters = true;}

		virtual double eval( std::vector<double>& args );
	};

//=================== Functions =====================================================

	inline double FickCoefficientProp::eval( std::vector<double>& args ){
		double result = d_params[d_nparams-1];
		for (size_t i= d_nparams-1; i>0; i--) {
			result = result * args[0] + d_params[i-1];
		}
		return result;
	}

}

//=================== Materials =====================================================

RadialPolynomial::RadialPolynomial()
{
		d_propertyMap = new std::map<std::string, boost::shared_ptr<AMP::Materials::Property<double> > >();
		INSERT_PROPERTY_IN_MAP( FickCoefficient, 				RadialPolynomial_NS);
}

}
}
