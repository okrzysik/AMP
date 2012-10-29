#ifndef included_AMP_DiffusionNonlinearFEOperatorParameters
#define included_AMP_DiffusionNonlinearFEOperatorParameters

#include <vector>

#include "vectors/Vector.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/diffusion/DiffusionConstants.h"


namespace AMP {
namespace Operator {

class DiffusionNonlinearFEOperatorParameters: public FEOperatorParameters {
public:

    DiffusionNonlinearFEOperatorParameters(const boost::shared_ptr<
            AMP::Database> &db) :
        FEOperatorParameters(db)
    {
    }

    virtual ~DiffusionNonlinearFEOperatorParameters() {
    }

    boost::shared_ptr<DiffusionTransportModel> d_transportModel;

    AMP::LinearAlgebra::Vector::shared_ptr d_FrozenTemperature;

    AMP::LinearAlgebra::Vector::shared_ptr d_FrozenConcentration;

    AMP::LinearAlgebra::Vector::shared_ptr d_FrozenBurnup;
protected:

private:

};

}
}

#endif

