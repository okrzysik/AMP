
#ifndef included_AMP_PelletStackMechanicsSolverParameters
#define included_AMP_PelletStackMechanicsSolverParameters

#include "solvers/ColumnSolver.h"

namespace AMP {
  namespace Solver {

    class PelletStackMechanicsSolverParameters : public SolverStrategyParameters {
      public:
        PelletStackMechanicsSolverParameters(){}
        PelletStackMechanicsSolverParameters(const boost::shared_ptr<AMP::Database> db) : SolverStrategyParameters(db) { }
        ~PelletStackMechanicsSolverParameters(){}

        boost::shared_ptr<AMP::Solver::ColumnSolver> d_columnSolver;

      protected:
      private:

    };

  }
}

#endif


