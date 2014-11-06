
#ifndef included_AMP_FirstOperator
#define included_AMP_FirstOperator

#include "operators/newFrozenVectorDesign/OnePointOperator.h"

namespace AMP {
  namespace Operator {

    class FirstOperator : public OnePointOperator {
      public :
        FirstOperator(const AMP::shared_ptr<OperatorParameters> & params) : OnePointOperator(params) {
          d_constant = 2.0;
          d_var.reset(new AMP::LinearAlgebra::Variable(params->d_db->getString("Variable")));
        }

        void apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u,
            AMP::LinearAlgebra::Vector::shared_ptr r, const double a = -1.0, const double b = 1.0) {
          AMP::LinearAlgebra::Vector::const_shared_ptr in = u->constSubsetVectorForVariable(d_var);
          AMP::LinearAlgebra::Vector::shared_ptr out = r->subsetVectorForVariable(d_var);
          out->scale((d_constant*a), in);
          if(f != NULL) {
            AMP::LinearAlgebra::Vector::const_shared_ptr rhs = f->constSubsetVectorForVariable(d_var);
            if(rhs != NULL) {
              out->axpy(b, rhs, out);
            }
          }
        }

        AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() {
          return d_var;
        }

        AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
          return d_var;
        }

      protected :
        AMP::LinearAlgebra::Variable::shared_ptr d_var;

      private :
    };

  }
}

#endif



