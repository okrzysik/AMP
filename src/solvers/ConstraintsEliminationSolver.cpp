#include "AMP/solvers/ConstraintsEliminationSolver.h"
#include "AMP/operators/ConstraintsEliminationOperator.h"

namespace AMP {
namespace Solver {

ConstraintsEliminationSolver::ConstraintsEliminationSolver(
    std::shared_ptr<ConstraintsEliminationSolverParameters> params )
    : SolverStrategy( params )
{
}

void ConstraintsEliminationSolver::apply( std::shared_ptr<const AMP::LinearAlgebra::Vector>,
                                          std::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    std::shared_ptr<AMP::Operator::ConstraintsEliminationOperator> op =
        std::dynamic_pointer_cast<AMP::Operator::ConstraintsEliminationOperator>( d_pOperator );
    op->copyMasterToSlave( u );
}

} // end namespace Solver
} // end namespace AMP
