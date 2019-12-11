
#ifndef included_AMP_MassLinearFEOperator
#define included_AMP_MassLinearFEOperator

// AMP files
#include "AMP/operators/libmesh/LinearFEOperator.h"
#include "AMP/operators/libmesh/MassLinearElement.h"
#include "AMP/operators/libmesh/MassLinearFEOperatorParameters.h"
#include "AMP/utils/Utilities.h"
#include <memory>

#include <vector>

namespace AMP {
namespace Operator {

class MassLinearFEOperator : public LinearFEOperator
{
public:
    explicit MassLinearFEOperator( const std::shared_ptr<MassLinearFEOperatorParameters> &params );

    virtual ~MassLinearFEOperator() {}

    void preAssembly( const std::shared_ptr<AMP::Operator::OperatorParameters> & ) override;

    void postAssembly() override;

    void preElementOperation( const AMP::Mesh::MeshElement & ) override;

    void postElementOperation() override;

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override;

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override;

    std::shared_ptr<MassDensityModel> getDensityModel() { return d_densityModel; };

protected:
    bool d_useConstantTemperature;

    bool d_useConstantConcentration;

    bool d_useConstantBurnup;

    double d_constantTemperatureValue;

    double d_constantConcentrationValue;

    double d_constantBurnupValue;

    AMP::LinearAlgebra::Vector::shared_ptr d_temperature;

    AMP::LinearAlgebra::Vector::shared_ptr d_concentration;

    AMP::LinearAlgebra::Vector::shared_ptr d_burnup;

    std::vector<std::vector<double>> d_elementMassMatrix;

    std::shared_ptr<MassLinearElement> d_massLinElem;

    std::shared_ptr<MassDensityModel> d_densityModel;

    std::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable;

    std::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;

private:
};
} // namespace Operator
} // namespace AMP

#endif
