#include "AMP/solvers/ConstraintsEliminationSolver.h"
#include "AMP/operators/ConstraintsEliminationOperator.h"

namespace AMP {
namespace Solver {

ConstraintsEliminationSolver::ConstraintsEliminationSolver(
    AMP::shared_ptr<ConstraintsEliminationSolverParameters> params )
    : SolverStrategy( params )
{
}

void ConstraintsEliminationSolver::solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector>,
                                          AMP::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    AMP::shared_ptr<AMP::Operator::ConstraintsEliminationOperator> op =
        AMP::dynamic_pointer_cast<AMP::Operator::ConstraintsEliminationOperator>( d_pOperator );
    op->copyMasterToSlave( u );
}

} // end namespace Solver
} // end namespace AMP
