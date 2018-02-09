#include "AMP/solvers/contact/MPCSolver.h"
#include "AMP/operators/contact/NodeToSegmentConstraintsOperator.h"

namespace AMP {
namespace Solver {

void MPCSolver::solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector>,
                       AMP::shared_ptr<AMP::LinearAlgebra::Vector> u )
{

    AMP::shared_ptr<AMP::Operator::NodeToSegmentConstraintsOperator> op =
        AMP::dynamic_pointer_cast<AMP::Operator::NodeToSegmentConstraintsOperator>( d_pOperator );

    op->applySolutionCorrection( u );
}
} // namespace Solver
} // namespace AMP
