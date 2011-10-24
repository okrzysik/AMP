#ifndef included_AMP_FickSoretNonlinearFEOperator
#define included_AMP_FickSoretNonlinearFEOperator

/* AMP files */
#include "boost/shared_ptr.hpp"
#include "ampmesh/MeshManager.h"
#include "ampmesh/DOFMap.h"
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

    typedef boost::shared_ptr <FickSoretNonlinearFEOperator>  shared_ptr;

    FickSoretNonlinearFEOperator(const boost::shared_ptr<OperatorParameters> & params);

    virtual ~FickSoretNonlinearFEOperator(){ }

    virtual void reset(const boost::shared_ptr<OperatorParameters>& params)
    {
        boost::shared_ptr<FickSoretNonlinearFEOperatorParameters> fsParams =
            boost::dynamic_pointer_cast<FickSoretNonlinearFEOperatorParameters>(params);

        d_FickOperator->reset(fsParams->d_FickParameters);
        d_SoretOperator->reset(fsParams->d_SoretParameters);
    }

    virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
            const  AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r,
            const double a = -1.0, const double b = 1.0);

    virtual boost::shared_ptr<OperatorParameters>
      getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& u)
    {
        return d_FickOperator->getJacobianParameters(u);
    }

    virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
      return d_OutputVariable;
    }

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1) {
        return d_FickOperator->getInputVariable(varId);
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

