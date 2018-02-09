#ifndef included_AMP_TrilinosThyraModelEvaluatorParameters
#define included_AMP_TrilinosThyraModelEvaluatorParameters


#include "AMP/operators/Operator.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/trilinos/nox/PrePostOperator.h"
#include "AMP/vectors/Vector.h"


namespace AMP {
namespace Solver {


/**
 * The TrilinosThyraModelEvaluator is a wrapper for a Thyra ModelEvaluator to
 * wrap AMP::Operators for use with Trilinos NOX solvers.
 */
class TrilinosThyraModelEvaluatorParameters
{
public:
    AMP::LinearAlgebra::Vector::shared_ptr d_icVec;             //!< The dofs to use for the vectors
    AMP::Operator::Operator::shared_ptr d_nonlinearOp;          //!< The non-linear operator
    AMP::Operator::Operator::shared_ptr d_linearOp;             //!< The linear operator
    AMP::Solver::SolverStrategy::shared_ptr d_preconditioner;   //!< The preconditioner
    AMP::Solver::PrePostOperator::shared_ptr d_prePostOperator; //!< The pre-post operator
};
} // namespace Solver
} // namespace AMP

#endif
