#ifndef included_AMP_TrilinosThyraModelEvaluatorParameters
#define included_AMP_TrilinosThyraModelEvaluatorParameters


#include "AMP/operators/Operator.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/trilinos/nox/PrePostOperator.h"
#include "AMP/vectors/Vector.h"


namespace AMP::Solver {


/**
 * The TrilinosThyraModelEvaluator is a wrapper for a Thyra ModelEvaluator to
 * wrap AMP::Operators for use with Trilinos NOX solvers.
 */
class TrilinosThyraModelEvaluatorParameters
{
public:
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_icVec;    //!< The dofs to use for the vectors
    std::shared_ptr<AMP::Operator::Operator> d_nonlinearOp; //!< The non-linear operator
    std::shared_ptr<AMP::Operator::Operator> d_linearOp;    //!< The linear operator
    std::shared_ptr<AMP::Solver::SolverStrategy> d_preconditioner;   //!< The preconditioner
    std::shared_ptr<AMP::Solver::PrePostOperator> d_prePostOperator; //!< The pre-post operator
};
} // namespace AMP::Solver

#endif
