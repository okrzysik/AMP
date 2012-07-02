/*
 * Property.cc
 *
 *  Created on: Sep 20, 2011
 *      Author: gad
 */

#include "Property.h"
#include "utils/Utilities.h"
#include "vectors/Vector.h"

#include <algorithm>

namespace AMP {
namespace Materials {

template<>
std::map<std::string, boost::shared_ptr<AMP::LinearAlgebra::Vector> >
Property<double>::make_map(const boost::shared_ptr<AMP::LinearAlgebra::MultiVector> &args)
{
	std::map<std::string, boost::shared_ptr<AMP::LinearAlgebra::Vector> > result;
	if (d_n_arguments>0) {
		size_t xls = d_translator.size();
		AMP_INSIST(xls>0, "attempt to make MultiVector map without setting translator");
		for (AMP::LinearAlgebra::MultiVector::vector_const_iterator vec=args->beginVector();
				vec != args->endVector(); ++vec)
		{
			std::string name = (*vec)->getVariable()->getName();

			for (std::map<std::string, std::string>::iterator pair=d_translator.begin();
					pair != d_translator.end(); pair++)
			{
				std::string key = pair->first;
				if (pair->second == name) {
					result.insert(std::make_pair(key, *vec));
				}
			}
		}
	}
	return result;
}

template<>
void Property<double>::evalv(boost::shared_ptr<AMP::LinearAlgebra::Vector>& r,
const std::map< std::string, boost::shared_ptr<AMP::LinearAlgebra::Vector> >& args)
{
	AMP_ASSERT(in_range(args));
	evalvActual(*r, args);
}

template<>
void Property<double>::evalv(boost::shared_ptr<AMP::LinearAlgebra::Vector>& r,
const boost::shared_ptr<AMP::LinearAlgebra::MultiVector>& args)
{
	std::map<std::string, boost::shared_ptr<AMP::LinearAlgebra::Vector> > mapargs = make_map(args);
	evalv(r, mapargs);
}

template<>
double Property<double>::NewtonSolve(double guess, double param1, double param2)
{
	double x_new = guess;
	double x_old = guess;
	bool converged = false;
	for (unsigned int iter=1; iter<=Newton_maxIter; ++iter){
		x_old = x_new;
		double perturbation = 1.0e-6;
		// numerical Jacobian with forward perturbation
		double J = (Residual(x_old+perturbation,param1,param2) - Residual(x_old,param1,param2))/perturbation;
		double dx = -1.0*Residual(x_old,param1,param2)/J;
		x_new = x_old + dx;
		// check convergence
		double abs_err = std::abs(x_new - x_old); // absolute error
		double rel_err = 0.0; // relative error
		if (x_old != 0.0){ // test to ensure no division by zero
			rel_err = std::abs((x_new - x_old)/x_old);
		}
		if ((abs_err < Newton_atol) and (rel_err < Newton_rtol)){
			converged = true;
			break;
		}
	}
	if (!converged){
		AMP_ERROR("Newton solve failed to converge for property function evaluation.");
	}
	return x_new;
}

template<>
double Property<double>::Residual(double arg1, double arg2, double arg3)
{
	AMP_ERROR("The function ``Property::Residual'' is virtual and expected to be overridden");
	return 0.0;
}

}
}
