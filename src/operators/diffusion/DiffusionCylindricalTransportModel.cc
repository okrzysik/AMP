/*
 * DiffusionCylindricalTransportModel.cc
 *
 *  Created on: Aug 19, 2011
 *      Author: gad
 */

#include "DiffusionCylindricalTransportModel.h"
#include "utils/Database.h"
#include <cmath>
#include "utils/Utilities.h"
#include <map>

/* Libmesh files */
#include "point.h"

namespace AMP {
namespace Operator {

DiffusionCylindricalTransportModel::
DiffusionCylindricalTransportModel(const boost::shared_ptr<DiffusionCylindricalTransportModelParameters> params):
DiffusionTransportTensorModel(params)
{
}

void
DiffusionCylindricalTransportModel::
getTensorTransport(
	std::vector< std::vector< std::vector<double> > >& result,
	std::map<std::string, boost::shared_ptr<std::vector<double> > >& args,
    const std::vector<Point>& coordinates)
{
	boost::shared_ptr<std::vector<double> >scaledp;
	double lower,upper;

	if (d_UseBilogScaling) {
		// do the transform
	    lower = d_BilogRange[0]+d_BilogEpsilonRangeLimit;
	    upper = d_BilogRange[1]-d_BilogEpsilonRangeLimit;
	    scaledp = bilogTransform(*args[d_BilogVariable], lower,upper);
		std::vector<double> &scaled = *scaledp;

		// save untransformed argument value
		for (size_t i=0; i<args[d_BilogVariable]->size(); i++) {
			double temp = (*args[d_BilogVariable])[i];
			(*args[d_BilogVariable])[i] = scaled[i];
			scaled[i] = temp;
		}
	}

	// evaluate material property as a function of radius
	// first fill in radius array
	std::vector<double> radialCoefficient(result.size());
	if (args.find("radius") != args.end()) {
		std::vector<double> &radius(*args["radius"]);
		for (size_t k=0; k<radius.size(); k++) {
			double x = coordinates[k](0);
			double y = coordinates[k](1);
			double r = sqrt(x*x+y*y);
			radius[k] = r;
		}
	}
	d_property->evalv(radialCoefficient, args);

	if (d_UseBilogScaling) {
		// restore untransformed argument value
		std::vector<double> &scaled = *scaledp;
		lower = d_BilogRange[0]+d_BilogEpsilonRangeLimit;
		upper = d_BilogRange[1]-d_BilogEpsilonRangeLimit;
		for (size_t i=0; i<(*args[d_BilogVariable]).size(); i++) {
			(*args[d_BilogVariable])[i] = scaled[i];
		}

		if (d_BilogScaleCoefficient) {
			for (size_t i=0; i<radialCoefficient.size(); i++) {
				bilogScale(radialCoefficient, lower,upper);
			}
		}
	}

	// form angle-dependent tensor factor
	AMP_ASSERT(coordinates.size() == result.size());
	for (size_t k=0; k<coordinates.size(); k++) {
		for (size_t i=0; i<3; i++) for (size_t j=0; j<3; j++) {
			result[k][i][j] = 0.;
		}
		double x = coordinates[k](0);
		double y = coordinates[k](1);
		double r2 = x*x+y*y;
		result[k][0][0] = x*x/r2;
		result[k][0][1] = x*y/r2;
		result[k][1][0]	= x*y/r2;
		result[k][1][1] = y*y/r2;
		for (size_t i=0; i<2; i++) for (size_t j=0; j<2; j++) {
			result[k][i][j] *= radialCoefficient[k];
		}
	}
}


}
}
