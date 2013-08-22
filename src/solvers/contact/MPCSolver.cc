
#include <solvers/contact/MPCSolver.h>
#include <operators/contact/NodeToSegmentConstraintsOperator.h>

namespace AMP {
  namespace Solver {

    void MPCSolver :: solve(boost::shared_ptr<const AMP::LinearAlgebra::Vector>, boost::shared_ptr<AMP::LinearAlgebra::Vector> u) {

      boost::shared_ptr<AMP::Operator::NodeToSegmentConstraintsOperator> op = boost::dynamic_pointer_cast<
        AMP::Operator::NodeToSegmentConstraintsOperator>(d_pOperator);

      op->applySolutionCorrection(u);

    }

  }
}



