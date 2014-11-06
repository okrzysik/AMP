
#ifndef included_contact_MPCSolver
#define included_contact_MPCSolver

#include <solvers/SolverStrategy.h>

namespace AMP {
  namespace Solver {

    typedef SolverStrategyParameters MPCSolverParameters;

    class MPCSolver : public SolverStrategy {
      public:
        MPCSolver(AMP::shared_ptr<MPCSolverParameters> params)
          : SolverStrategy(params) { }

        void solve(AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f, AMP::shared_ptr<AMP::LinearAlgebra::Vector> u);

    };

  }
}

#endif


