
#ifndef included_AMP_ConstraintsEliminationSolver
#define included_AMP_ConstraintsEliminationSolver

#include <solvers/SolverStrategy.h>

namespace AMP {
  namespace Solver {

    typedef SolverStrategyParameters ConstraintsEliminationSolverParameters;

    class ConstraintsEliminationSolver : public SolverStrategy {
      public:
        explicit ConstraintsEliminationSolver(AMP::shared_ptr<ConstraintsEliminationSolverParameters> params);
        void solve(AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f, AMP::shared_ptr<AMP::LinearAlgebra::Vector> u);
    };

  }
}

#endif


