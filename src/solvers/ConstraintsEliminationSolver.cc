
#include <solvers/ConstraintsEliminationSolver.h>
#include <operators/ConstraintsEliminationOperator.h>

namespace AMP {
  namespace Solver {
    
    ConstraintsEliminationSolver::ConstraintsEliminationSolver(boost::shared_ptr<ConstraintsEliminationSolverParameters> params)
      : SolverStrategy(params) 
    { 

    }

    void ConstraintsEliminationSolver::solve(boost::shared_ptr<AMP::LinearAlgebra::Vector>, 
      boost::shared_ptr<AMP::LinearAlgebra::Vector> u) 
    {
      boost::shared_ptr<AMP::Operator::ConstraintsEliminationOperator> op = boost::dynamic_pointer_cast<
        AMP::Operator::ConstraintsEliminationOperator>(d_pOperator);
      op->copyMasterToSlave(u);
    }

  } // end namespace Solver
} // end namespace AMP

