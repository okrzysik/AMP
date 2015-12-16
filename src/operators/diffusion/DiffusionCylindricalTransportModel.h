/*
 * DiffusionCylindricalTransportModel.h
 *
 *  Created on: Aug 19, 2011
 *      Author: gad
 */

#ifndef DIFFUSIONCYLINDRICALTRANSPORTMODEL_H_
#define DIFFUSIONCYLINDRICALTRANSPORTMODEL_H_

#include "operators/diffusion/DiffusionTransportTensorModel.h"
#include "utils/Utilities.h"
#include <string>

namespace AMP {
namespace Operator {
typedef ElementPhysicsModelParameters DiffusionCylindricalTransportModelParameters;

class DiffusionCylindricalTransportModel : public DiffusionTransportTensorModel
{
public:
    explicit DiffusionCylindricalTransportModel(
        const AMP::shared_ptr<DiffusionTransportTensorModelParameters> params );

    /**
     * \brief transport model returning a vector of tensors for cylindrical symmetry
     * \param result result[i] is a tensor of diffusion coefficients.
     * \param args args[j][i] is j-th material evalv argument
     * \param Coordinates vector of points that define the spatial location
     */
    virtual void
    getTensorTransport( std::vector<std::vector<AMP::shared_ptr<std::vector<double>>>> &result,
                        std::map<std::string, AMP::shared_ptr<std::vector<double>>> &args,
                        const std::vector<libMesh::Point> &Coordinates = d_DummyCoords );

private:
    std::string d_RadiusArgument;
};
}
}

#endif /* DIFFUSIONCYLINDRICALTRANSPORTMODEL_H_ */
