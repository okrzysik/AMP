/*
 * CylindricallySymmetric.h
 *
 *  Created on: June 11, 2010
 *	  Author: bm
 */

#include "CylindricallySymmetric.h"

#include "materials/Property.h"
#include "materials/VectorProperty.h"
#include "materials/TensorProperty.h"
#include "materials/Material.h"
#include "utils/Utilities.h"

#include <string>
#include <valarray>
#include <cmath>
#include <limits>
#include <iostream>


static inline int lround(double x) { return x>=0.0 ? ceil(x):floor(x); }

namespace AMP
{
namespace Materials
{

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
namespace CylindricallySymmetric_NS
{

//=================== Constants =====================================================

static const std::string name_base("CylindricallySymmetric");
static const std::string source("");

double Pi=3.1415926535898;

static const double rMinVal =  0.0;
static const double rMaxVal =  std::numeric_limits<double>::max();
static const double tMinVal =  0.0;
static const double tMaxVal =  2.*Pi;
static const double zMinVal = -std::numeric_limits<double>::max();
static const double zMaxVal =  std::numeric_limits<double>::max();

static const std::string argumentsRadialFick[1] =
{ "radius" };
static const unsigned int nargumentsRadialFick = 1;
static const double fickRadialRanges[1][2] =
{
{ rMinVal, rMaxVal } };

static const std::string argumentsLongitudinalFick[1] =
{ "zee" };
static const unsigned int nargumentsLongitudinalFick = 1;
static const double fickLongitudinalRanges[1][2] =
{
{ zMinVal, zMaxVal } };

static const std::string argumentsTensorFick[3] =
{ "radius", "theta", "zee" };
static const unsigned int nargumentsTensorFick = 3;
static const double fickTensorRanges[3][2] =
{
{ rMinVal, rMaxVal },
{ tMinVal, tMaxVal },
{ zMinVal, zMaxVal } };

// default polynomial is a constant
static const size_t nparams = 1;
static const double params[] =
{ 1. };
static const size_t nparamsTensor = 3;
static const double paramsTensor[] =
{ 1., 1., 1. };

//=================== Classes =======================================================

/** radial diffusion coefficient */
class ScalarRadialFickProp: public Property<double>
{
public:
	ScalarRadialFickProp() :
			Property<double>(name_base + "_" + "ScalarRadialFick", // Name string
			source, // Reference source
					params, // Property parameters
					nparams, // Number of parameters
					argumentsRadialFick, // Names of arguments to the eval function
					nargumentsRadialFick, // Number of arguments
					fickRadialRanges) // Ranges
	{
		d_variableNumberParameters = true;
	}

	virtual void set_parameters_and_number(const double* params,
			const unsigned int nparams)
	{
		AMP_ASSERT(nparams>0);
		Property<double>::set_parameters_and_number(params, nparams);
	}

	/** returns property and derivative wrto r
	 * \return [0]=property, [1]=property derivative wrto r
	 */
	virtual double eval(std::vector<double>& args);
};

/** radial diffusion coefficient */
class RadialFickProp: public VectorProperty<double>
{
public:
	RadialFickProp() :
			VectorProperty<double>(name_base + "_" + "RadialFick", // Name string
			source, // Reference source
					params, // Property parameters
					nparams, // Number of parameters
					argumentsRadialFick, // Names of arguments to the eval function
					nargumentsRadialFick, // Number of arguments
					fickRadialRanges) // Ranges
	{
		d_variableNumberParameters = true;
	}

	virtual void set_parameters_and_number(const double* params,
			const unsigned int nparams)
	{
		AMP_ASSERT(nparams>0);
		Property<double>::set_parameters_and_number(params, nparams);
	}

	/** returns property and derivative wrto r
	 * \return [0]=property, [1]=property derivative wrto r
	 */
	virtual std::vector<double>
	evalVector(std::vector<double>& args);
};

/** longitudinal diffusion coefficient */
class LongitudinalFickProp: public VectorProperty<double>
{
public:
	LongitudinalFickProp() :
			VectorProperty<double>(name_base + "_" + "LongitudinalFick", // Name string
			source, // Reference source
					params, // Property parameters
					nparams, // Number of parameters
					argumentsLongitudinalFick, // Names of arguments to the eval function
					nargumentsLongitudinalFick, // Number of arguments
					fickLongitudinalRanges) // Ranges
	{
		d_variableNumberParameters = true;
	}

	virtual void set_parameters_and_number(const double* params,
			const unsigned int nparams)
	{
		AMP_ASSERT(nparams>0);
		Property<double>::set_parameters_and_number(params, nparams);
	}

	virtual std::vector<double>
	evalVector(std::vector<double>& args);
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
class TensorFickProp: public TensorProperty<double>
{
public:
	TensorFickProp() :
			TensorProperty<double>(name_base + "_" + "TensorFick", // Name string
			source, // Reference source
					paramsTensor, // Property parameters
					nparamsTensor, // Number of parameters
					argumentsTensorFick, // Names of arguments to the eval function
					nargumentsTensorFick, // Number of arguments
					fickTensorRanges, // ranges
					std::vector<size_t>(2, 3)) // dimensions
	{
		d_variableNumberParameters = true;
		d_variableDimensions = true;
		d_AuxiliaryDataInteger.insert(std::make_pair("derivative", 0));
		set_parameters_and_number(paramsTensor, nparamsTensor);
	}

	// NOTE: must change dimension first before changing number of parameters
	virtual void set_parameters_and_number(const double* params,
			const unsigned int nparams)
	{
		Property<double>::set_parameters_and_number(params, nparams);
		AMP_ASSERT(d_nparams >= 3);
		d_nparamsRadial = lround(d_params[0]);
		AMP_ASSERT(d_nparamsRadial < d_nparams-1);
		d_nparamsLongitudinal = d_nparams - 1 - d_nparamsRadial;
		d_radialK.set_parameters_and_number(&d_params[1], d_nparamsRadial);
		d_longitudinalK.set_parameters_and_number(&d_params[1+d_nparamsRadial], d_nparamsLongitudinal);
	}

	virtual std::vector<std::vector<double> >
	evalTensor(std::vector<double>& args);

private:
	RadialFickProp d_radialK;
	LongitudinalFickProp d_longitudinalK;
	unsigned int d_nparamsRadial;
	unsigned int d_nparamsLongitudinal;
};

//=================== Functions =====================================================

inline double ScalarRadialFickProp::eval(std::vector<double>& args)
{
	AMP_ASSERT(!args.empty());
	double result;
	result = d_params[d_nparams - 1];
	for (size_t i = d_nparams - 1; i > 0; i--)
	{
		result = result * args[0] + d_params[i - 1];
	}
	return result;
}

inline std::vector<double> RadialFickProp::evalVector(std::vector<double>& args)
{
	AMP_ASSERT(!args.empty());
	std::vector<double> result(2);
	result[0] = d_params[d_nparams - 1];
	result[1] = (d_nparams-1)*d_params[d_nparams - 1];
	for (size_t i = d_nparams - 1; i > 0; i--)
	{
		result[0] = result[0] * args[0] + d_params[i - 1];
	}
	for (size_t i = d_nparams - 1; i > 1; i--)
	{
		result[1] = result[1] * args[0] + (i-1)*d_params[i - 1];
	}
	return result;
}

inline std::vector<double> LongitudinalFickProp::evalVector(std::vector<double>& args)
{
	AMP_ASSERT(!args.empty());
	std::vector<double> result(2);
	result[0] = d_params[d_nparams - 1];
	result[1] = (d_nparams-1)*d_params[d_nparams - 1];
	for (size_t i = d_nparams - 1; i > 0; i--)
	{
		result[0] = result[0] * args[0] + d_params[i - 1];
	}
	for (size_t i = d_nparams - 1; i > 1; i--)
	{
		result[1] = result[1] * args[0] + (i-1)*d_params[i - 1];
	}
	return result;
}

std::vector<std::vector<double> > TensorFickProp::evalTensor(
		std::vector<double>& args)
{
	AMP_ASSERT(args.size()>2);
	std::vector<std::vector<double> > result(3,std::vector<double>(3,0.));
	std::vector<double> argr(1,args[0]);
	std::vector<double> argz(1,args[2]);
	std::vector<double> Kr = d_radialK.evalVector(argr);
	std::vector<double> Kz = d_longitudinalK.evalVector(argz);
	double cth = cos(args[1]);
	double sth = sin(args[1]);
	int deriv = d_AuxiliaryDataInteger.find("derivative")->second;
	switch(deriv) {
	case 0:
	{
		result[0][0] = cth*cth*Kr[0];
		result[0][1] = sth*cth*Kr[0];
		result[1][0] = sth*cth*Kr[0];
		result[1][1] = sth*sth*Kr[0];
		result[2][2] =         Kz[0];
		break;
	}
	case 1:
	{
		result[0][0] = cth*cth*Kr[1];
		result[0][1] = sth*cth*Kr[1];
		result[1][0] = sth*cth*Kr[1];
		result[1][1] = sth*sth*Kr[1];
		result[2][2] = 0.;
		break;
	}
	case 2:
	{
		result[0][0] = 0.;
		result[0][1] = 0.;
		result[1][0] = 0.;
		result[1][1] = 0.;
		result[2][2] = Kz[1];
		break;
	}
	default:
		AMP_ASSERT(false);
		break;
	}
	return result;
}

}

//=================== Materials =====================================================

CylindricallySymmetric::CylindricallySymmetric()
{
	d_propertyMap = new std::map<std::string,
			boost::shared_ptr<AMP::Materials::Property<double> > >();
	INSERT_PROPERTY_IN_MAP( ScalarRadialFick, CylindricallySymmetric_NS);
	INSERT_PROPERTY_IN_MAP( RadialFick,       CylindricallySymmetric_NS);
	INSERT_PROPERTY_IN_MAP( LongitudinalFick, CylindricallySymmetric_NS);
	INSERT_PROPERTY_IN_MAP( TensorFick,       CylindricallySymmetric_NS);
}

}
}
