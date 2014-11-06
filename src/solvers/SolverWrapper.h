
#ifndef included_SolverWrapper
#define included_SolverWrapper

#include "solvers/SolverStrategy.h"
#include "operators/OperatorWrapper.h"

namespace AMP {
  namespace Solver {

    class SolverWrapper : public SolverStrategy {
      public:
        SolverWrapper() { }

        virtual ~SolverWrapper() { }

        void setSolver(AMP::shared_ptr<SolverStrategy> in) {
          d_solver = in;
          d_pOperator = in->getOperator();
        }

        void solve(AMP::shared_ptr<const AMP::LinearAlgebra::Vector>  f,
            AMP::shared_ptr<AMP::LinearAlgebra::Vector>  u) {
          AMP::shared_ptr<AMP::Operator::OperatorWrapper> wrapOp = AMP::dynamic_pointer_cast<
            AMP::Operator::OperatorWrapper>(d_pOperator);
          wrapOp->setFullVector(u);
          AMP::Vector::shared_ptr uTmp = u->subsetVectorForVariable(d_pOperator->getOutputVariable());
          d_solver->solve(f, uTmp);
        }

      protected:
        AMP::shared_ptr<SolverStrategy> d_solver;
    };

  }
}

#endif


