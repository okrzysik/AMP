
#include "LinearOperator.h"
#include "utils/Utilities.h"

namespace AMP {
  namespace Operator {

    LinearOperator :: LinearOperator (const boost::shared_ptr<OperatorParameters> & params)
      : Operator (params) {
        d_matrix.reset();
        d_applyCount = 0;
      }

    boost::shared_ptr<AMP::LinearAlgebra::Matrix> LinearOperator :: getMatrix() {
      return d_matrix;
    }

    void LinearOperator :: setMatrix(const boost::shared_ptr<AMP::LinearAlgebra::Matrix> & in_mat) {
      d_matrix = in_mat;
    }

    void LinearOperator :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
        AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a, const double b)
    {
      d_applyCount++;

      AMP_INSIST( ((u.get()) != NULL), "NULL Solution Vector" );
      AMP_INSIST( ((r.get()) != NULL), "NULL Residual Vector" );
      AMP_INSIST( ((d_matrix.get()) != NULL), "NULL Matrix" );

      AMP::LinearAlgebra::Variable::shared_ptr tmpInpVar = getInputVariable();
      AMP::LinearAlgebra::Variable::shared_ptr tmpOutVar = getOutputVariable();

      AMP::LinearAlgebra::Vector::shared_ptr uInternal = u->subsetVectorForVariable( tmpInpVar );
      AMP::LinearAlgebra::Vector::shared_ptr rInternal = r->subsetVectorForVariable( tmpOutVar );

      AMP_INSIST( (uInternal.get() != NULL), "uInternal is NULL" );
      AMP_INSIST( (rInternal.get() != NULL), "rInternal is NULL" );

      d_matrix->mult(uInternal, rInternal);

      if(f.get() == NULL) {
        rInternal->scale(a);
      } else {
        AMP::LinearAlgebra::Vector::shared_ptr fInternal = f->subsetVectorForVariable( tmpOutVar );
        if(fInternal.get() == NULL) {
          rInternal->scale(a);
        } else {
          rInternal->axpby(b, a, fInternal);
        }
      }
    }

  }
}

