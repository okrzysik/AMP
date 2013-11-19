#ifndef included_AMP_TrilinosThyraModelEvaluatorParameters
#define included_AMP_TrilinosThyraModelEvaluatorParameters


#include "vectors/Vector.h"
#include "operators/Operator.h"
#include "solvers/SolverStrategy.h"


namespace AMP {
namespace Solver {


/**
  * The TrilinosThyraModelEvaluator is a wrapper for a Thyra ModelEvaluator to 
  * wrap AMP::Operators for use with Trilinos NOX solvers.
  */
class TrilinosThyraModelEvaluatorParameters
{
public:
    
    AMP::LinearAlgebra::Vector::shared_ptr  d_icVec;            //!< The dofs to use for the vectors
    AMP::Operator::Operator::shared_ptr     d_nonlinearOp;      //!< The non-linear operator
    AMP::Operator::Operator::shared_ptr     d_linearOp;         //!< The linear operator
    AMP::Solver::SolverStrategy::shared_ptr d_preconditioner;   //!< The preconditioner

};


}
}

#endif

