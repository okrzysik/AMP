#include "AMP/solvers/contact/MPCSolver.h"
#include "AMP/operators/contact/NodeToSegmentConstraintsOperator.h"

namespace AMP {
namespace Solver {

void MPCSolver::apply( std::shared_ptr<const AMP::LinearAlgebra::Vector>,
                       std::shared_ptr<AMP::LinearAlgebra::Vector> u )
{

    std::shared_ptr<AMP::Operator::NodeToSegmentConstraintsOperator> op =
        std::dynamic_pointer_cast<AMP::Operator::NodeToSegmentConstraintsOperator>( d_pOperator );

    op->applySolutionCorrection( u );
}
} // namespace Solver
} // namespace AMP
