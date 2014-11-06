#ifndef included_AMP_FickSoretNonlinearFEOperator
#define included_AMP_FickSoretNonlinearFEOperator

/* AMP files */
#include "utils/shared_ptr.h"
#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "discretization/DOF_Manager.h"
#include "operators/diffusion/FickSoretNonlinearFEOperatorParameters.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"

namespace AMP {
namespace Operator {

  /**
    Class to add the output of the Fick and Soret operators.
    */
  class FickSoretNonlinearFEOperator : public Operator
  {
  public :

    typedef AMP::shared_ptr <FickSoretNonlinearFEOperator>  shared_ptr;

    FickSoretNonlinearFEOperator(const AMP::shared_ptr<OperatorParameters> & params);

    virtual ~FickSoretNonlinearFEOperator(){ }

    virtual void reset(const AMP::shared_ptr<OperatorParameters>& params)
    {
        AMP::shared_ptr<FickSoretNonlinearFEOperatorParameters> fsParams =
            AMP::dynamic_pointer_cast<FickSoretNonlinearFEOperatorParameters>(params);

        d_FickOperator->reset(fsParams->d_FickParameters);
        d_SoretOperator->reset(fsParams->d_SoretParameters);
    }

    virtual void apply(AMP::LinearAlgebra::Vector::const_shared_ptr f,
            AMP::LinearAlgebra::Vector::const_shared_ptr u, AMP::LinearAlgebra::Vector::shared_ptr r,
            const double a = -1.0, const double b = 1.0);

    virtual AMP::shared_ptr<OperatorParameters>
      getJacobianParameters(const AMP::shared_ptr<AMP::LinearAlgebra::Vector>& u)
    {
        return d_FickOperator->getJacobianParameters(u);
    }

    virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
      return d_OutputVariable;
    }

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() {
        return d_FickOperator->getInputVariable();
    }

    /**
     * checks input to apply operator for satisfaction of range conditions
     */
    bool isValidInput(AMP::LinearAlgebra::Vector::shared_ptr &u){
        bool result;
        result = d_FickOperator->isValidInput(u) and d_SoretOperator->isValidInput(u);
        return result;
    }

    DiffusionNonlinearFEOperator::shared_ptr getFickOperator(){return d_FickOperator;}
    DiffusionNonlinearFEOperator::shared_ptr getSoretOperator(){return d_SoretOperator;}

  private:

    DiffusionNonlinearFEOperator::shared_ptr d_FickOperator;
    DiffusionNonlinearFEOperator::shared_ptr d_SoretOperator;

    AMP::LinearAlgebra::Variable::shared_ptr d_OutputVariable;

    bool d_AddSoretTerm;
  };

}
}

#endif

