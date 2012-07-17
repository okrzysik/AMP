
#ifndef included_MPCSolver
#define included_MPCSolver

#include "SolverStrategy.h"

namespace AMP {
  namespace Solver {

    typedef SolverStrategyParameters MPCSolverParameters;

    class MPCSolver : public SolverStrategy {
      public:
        MPCSolver(boost::shared_ptr<MPCSolverParameters> params) 
          : SolverStrategy(params) { }

        void solve(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, boost::shared_ptr<AMP::LinearAlgebra::Vector> u);

    };

  }
}

#endif


