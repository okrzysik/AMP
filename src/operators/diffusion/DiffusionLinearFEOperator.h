
#ifndef included_AMP_DiffusionLinearFEOperator
#define included_AMP_DiffusionLinearFEOperator

/* AMP files */
#include "AMP/ampmesh/MeshElement.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperatorParameters.h"
#include "AMP/operators/libmesh/LinearFEOperator.h"
#include "AMP/utils/Utilities.h"

/* Boost files */
#include "AMP/utils/shared_ptr.h"

#include <vector>


namespace AMP {
namespace Operator {


class DiffusionLinearFEOperator : public LinearFEOperator
{
public:
    explicit DiffusionLinearFEOperator(
        const AMP::shared_ptr<DiffusionLinearFEOperatorParameters> &params );

    virtual ~DiffusionLinearFEOperator() {}

    void preAssembly( const AMP::shared_ptr<OperatorParameters> &params ) override;

    void postAssembly() override;

    void preElementOperation( const AMP::Mesh::MeshElement & ) override;

    void postElementOperation() override;

    AMP::shared_ptr<DiffusionTransportModel> getTransportModel();

protected:
    bool d_useConstantTemperature;

    bool d_useConstantConcentration;

    bool d_useConstantBurnup;

    AMP::LinearAlgebra::Vector::shared_ptr d_temperature;

    AMP::LinearAlgebra::Vector::shared_ptr d_concentration;

    AMP::LinearAlgebra::Vector::shared_ptr d_burnup;

    std::vector<std::vector<double>> d_elementStiffnessMatrix;

    AMP::shared_ptr<DiffusionLinearElement> d_diffLinElem;

    AMP::shared_ptr<DiffusionTransportModel> d_transportModel;
};
} // namespace Operator
} // namespace AMP

#endif
