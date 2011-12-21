
#ifndef included_AMP_OnePointSolver
#define included_AMP_OnePointSolver

#include "solvers/SolverStrategy.h"
#include "operators/newFrozenVectorDesign/OnePointOperator.h"

namespace AMP {
  namespace Solver {

    class OnePointSolver: public SolverStrategy {
      public:
        OnePointSolver(boost::shared_ptr<SolverStrategyParameters> params) : SolverStrategy(params) {
          d_onePointOp = boost::dynamic_pointer_cast<AMP::Operator::OnePointOperator>(d_pOperator);
        }

        void solve(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, boost::shared_ptr<AMP::LinearAlgebra::Vector> u) {

        }

      protected:
        boost::shared_ptr<AMP::Operator::OnePointOperator> d_onePointOp;
    };

  }
}

#endif



