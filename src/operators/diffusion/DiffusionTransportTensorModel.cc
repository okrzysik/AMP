/*
 * DiffusionTransportModel.cc
 *
 *  Created on: Jul 8, 2010
 *      Author: gad
 */

#include "DiffusionTransportTensorModel.h"
#include "utils/Database.h"
#include <cmath>

namespace AMP {
namespace Operator {

DiffusionTransportTensorModel::
DiffusionTransportTensorModel(const boost::shared_ptr<DiffusionTransportTensorModelParameters> params):
DiffusionTransportModel(params)
{
	d_IsTensor = true;
}

void
DiffusionTransportTensorModel::
getTensorTransport(
    std::vector< std::vector< boost::shared_ptr<std::vector<double> > > >& result,
    std::map<std::string, boost::shared_ptr<std::vector<double> > >& args,
    const std::vector<Point>& Coordinates)
{
	boost::shared_ptr<std::vector<double> >scaledp;
	double lower,upper;

	if (d_UseBilogScaling) {
		// do the transform
	  lower = d_BilogRange[0]+d_BilogEpsilonRangeLimit;
	  upper = d_BilogRange[1]-d_BilogEpsilonRangeLimit;
	  scaledp = bilogTransform((*args[d_BilogVariable]), lower,upper);
		std::vector<double> &scaled = *scaledp;

		// save untransformed argument value
		for (size_t i=0; i<(*args[d_BilogVariable]).size(); i++) {
			double temp = (*args[d_BilogVariable])[i];
			(*args[d_BilogVariable])[i] = scaled[i];
			scaled[i] = temp;
		}
	}

	// evaluate material property
	// material library has been temporarily supplied with a dummy evalv for tensors
	// new material interface will fix.
	d_property->evalv(result, args);

	if (d_UseBilogScaling) {
		// restore untransformed argument value
		std::vector<double> &scaled = *scaledp;
		lower = d_BilogRange[0]+d_BilogEpsilonRangeLimit;
		upper = d_BilogRange[1]-d_BilogEpsilonRangeLimit;
		for (size_t i=0; i<(*args[d_BilogVariable]).size(); i++) {
			(*args[d_BilogVariable])[i] = scaled[i];
		}

		if (d_BilogScaleCoefficient) {
			for (size_t i=0; i<result.size(); i++) for (size_t j=0; j<result[i].size(); j++) {
				bilogScale(*result[i][j], lower,upper);
			}
		}
	}
}

}
}

