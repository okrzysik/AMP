#ifndef included_AMP_FickSoretNonlinearFEOperator
#define included_AMP_FickSoretNonlinearFEOperator

/* AMP files */
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshElement.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/operators/diffusion/FickSoretNonlinearFEOperatorParameters.h"
#include <memory>

namespace AMP {
namespace Operator {

/**
  Class to add the output of the Fick and Soret operators.
  */
class FickSoretNonlinearFEOperator : public Operator
{
public:
    typedef std::shared_ptr<FickSoretNonlinearFEOperator> shared_ptr;

    explicit FickSoretNonlinearFEOperator( const std::shared_ptr<OperatorParameters> &params );

    virtual ~FickSoretNonlinearFEOperator() {}

    std::string type() const override { return "FickSoretNonlinearFEOperator"; }

    void reset( const std::shared_ptr<OperatorParameters> &params ) override
    {
        std::shared_ptr<FickSoretNonlinearFEOperatorParameters> fsParams =
            std::dynamic_pointer_cast<FickSoretNonlinearFEOperatorParameters>( params );

        d_FickOperator->reset( fsParams->d_FickParameters );
        d_SoretOperator->reset( fsParams->d_SoretParameters );
    }

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    std::shared_ptr<OperatorParameters>
    getParameters( const std::string &type,
                   AMP::LinearAlgebra::Vector::const_shared_ptr u,
                   std::shared_ptr<OperatorParameters> params = nullptr ) override
    {
        return d_FickOperator->getParameters( type, u, params );
    }

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override
    {
        return d_OutputVariable;
    }

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override
    {
        return d_FickOperator->getInputVariable();
    }

    /**
     * checks input to apply operator for satisfaction of range conditions
     */
    bool isValidInput( AMP::LinearAlgebra::Vector::shared_ptr &u ) override
    {
        bool result;
        result = d_FickOperator->isValidInput( u ) and d_SoretOperator->isValidInput( u );
        return result;
    }

    DiffusionNonlinearFEOperator::shared_ptr getFickOperator() { return d_FickOperator; }
    DiffusionNonlinearFEOperator::shared_ptr getSoretOperator() { return d_SoretOperator; }

protected:
private:
    DiffusionNonlinearFEOperator::shared_ptr d_FickOperator;
    DiffusionNonlinearFEOperator::shared_ptr d_SoretOperator;

    AMP::LinearAlgebra::Variable::shared_ptr d_OutputVariable;

    bool d_AddSoretTerm;
};
} // namespace Operator
} // namespace AMP

#endif
