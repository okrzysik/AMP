
#include "operators/MoveMeshOperator.h"

namespace AMP {
  namespace Operator {

    MoveMeshOperator :: MoveMeshOperator(const boost::shared_ptr<OperatorParameters>& params)
      : Operator(params) {
        d_prevDisp.reset();
        d_var.reset();
      }

    void MoveMeshOperator :: setVariable(AMP::LinearAlgebra::Variable::shared_ptr var) {
      d_var = var;
    }

    AMP::LinearAlgebra::Variable::shared_ptr MoveMeshOperator :: getInputVariable(int varId) {
      return d_var;
    }

    void MoveMeshOperator :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, 
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r,
        const double a, const double b)  {
      AMP::LinearAlgebra::Vector::shared_ptr dispVec = u->subsetVectorForVariable(d_var);

      if(d_prevDisp == NULL) {
        d_prevDisp = dispVec->cloneVector();
        d_prevDisp->zero();
      }

      AMP::LinearAlgebra::Vector::shared_ptr deltaDisp = dispVec->cloneVector();
      deltaDisp->subtract(dispVec, d_prevDisp);

      d_Mesh->displaceMesh(deltaDisp);

      d_prevDisp->copyVector(dispVec);
    }

  }
}




