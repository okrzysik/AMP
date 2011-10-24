#ifndef included_AMP_DiffusionLinearFEOperatorParameters
#define included_AMP_DiffusionLinearFEOperatorParameters

#include "operators/FEOperatorParameters.h"

#include "operators/diffusion/DiffusionTransportModel.h"

namespace AMP {
namespace Operator {

class DiffusionLinearFEOperatorParameters: public FEOperatorParameters {
public:

    DiffusionLinearFEOperatorParameters(
            const boost::shared_ptr<AMP::Database> &db) :
        FEOperatorParameters(db) {
    }

    virtual ~DiffusionLinearFEOperatorParameters() {
    }

    boost::shared_ptr<DiffusionTransportModel> d_transportModel;

    AMP::LinearAlgebra::Vector::shared_ptr d_temperature;

    AMP::LinearAlgebra::Vector::shared_ptr d_concentration;

    AMP::LinearAlgebra::Vector::shared_ptr d_burnup;

protected:

private:

};

}
}

#endif

