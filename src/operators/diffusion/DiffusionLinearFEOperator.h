#ifndef included_AMP_DiffusionLinearFEOperator
#define included_AMP_DiffusionLinearFEOperator

#include "AMP/mesh/MeshElement.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperatorParameters.h"
#include "AMP/operators/libmesh/LinearFEOperator.h"
#include "AMP/utils/Utilities.h"
#include <memory>

#include <vector>


namespace AMP::Operator {


class DiffusionLinearFEOperator : public LinearFEOperator
{
public:
    explicit DiffusionLinearFEOperator( std::shared_ptr<const OperatorParameters> params );

    virtual ~DiffusionLinearFEOperator() {}

    void preAssembly( std::shared_ptr<const OperatorParameters> params ) override;

    void postAssembly() override;

    void preElementOperation( const AMP::Mesh::MeshElement & ) override;

    void postElementOperation() override;

    std::shared_ptr<DiffusionTransportModel> getTransportModel();

    void set_temperature( std::shared_ptr<const AMP::LinearAlgebra::Vector> );
    void set_concentration( std::shared_ptr<const AMP::LinearAlgebra::Vector> );
    void set_burnup( std::shared_ptr<const AMP::LinearAlgebra::Vector> );

protected:
    bool d_useConstantTemperature;
    bool d_useConstantConcentration;
    bool d_useConstantBurnup;

    std::shared_ptr<const AMP::LinearAlgebra::Vector> d_temperature;
    std::shared_ptr<const AMP::LinearAlgebra::Vector> d_concentration;
    std::shared_ptr<const AMP::LinearAlgebra::Vector> d_burnup;

    std::vector<std::vector<double>> d_elementStiffnessMatrix;

    std::shared_ptr<DiffusionLinearElement> d_diffLinElem;

    std::shared_ptr<DiffusionTransportModel> d_transportModel;
};
} // namespace AMP::Operator

#endif
