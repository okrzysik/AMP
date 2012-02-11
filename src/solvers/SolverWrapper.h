
#ifndef included_SolverWrapper
#define included_SolverWrapper

#include "solvers/SolverStrategy.h"
#include "operators/OperatorWrapper.h"

namespace AMP {
  namespace Solver {

    class SolverWrapper : public SolverStrategy {
      public:
        SolverWrapper(boost::shared_ptr<SolverStrategyParameters> params) 
          : SolverStrategy(params) { }

        virtual ~SolverWrapper() { }

        void setSolver(boost::shared_ptr<SolverStrategy> in) {
          d_solver = in;
        }

        void solve(boost::shared_ptr<AMP::LinearAlgebra::Vector>  f,
            boost::shared_ptr<AMP::LinearAlgebra::Vector>  u) {
          boost::shared_ptr<AMP::Operator::OperatorWrapper> wrapOp = boost::dynamic_pointer_cast<
            AMP::Operator::OperatorWrapper>(d_pOperator);
          wrapOp->setFullVector(u);
          AMP::Vector::shared_ptr uTmp = u->subsetVectorForVariable(d_pOperator->getOutputVariable());
          d_solver->solve(f, uTmp);
        }

      protected:
        boost::shared_ptr<SolverStrategy> d_solver;
    };

  }
}

#endif


