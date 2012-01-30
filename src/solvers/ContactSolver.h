
#ifndef included_ContactSolver
#define included_ContactSolver

#include "SolverStrategy.h"

namespace AMP {
  namespace Solver {

    typedef SolverStrategyParameters ContactSolverParameters;

    class ContactSolver : public SolverStrategy {
      public:
        ContactSolver() { }

        ContactSolver(boost::shared_ptr<ContactSolverParameters> params) : SolverStrategy(params) { }

        ~ContactSolver() { }

        void solve(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, boost::shared_ptr<AMP::LinearAlgebra::Vector>  u);

        void initialize(boost::shared_ptr<SolverStrategyParameters> const parameters) { }

        void setInitialGuess( boost::shared_ptr<AMP::LinearAlgebra::Vector>  initialGuess ) { } 

      private:

    };

  }
}

#endif


