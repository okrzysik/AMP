
#ifndef included_AMP_OnePointSolver
#define included_AMP_OnePointSolver

#include "operators/newFrozenVectorDesign/OnePointOperator.h"
#include "solvers/SolverStrategy.h"
#include "vectors/newFrozenVectorDesign/newFrozenVectorDesignHelpers.h"

namespace AMP {
namespace Solver {

class OnePointSolver : public SolverStrategy
{
public:
    explicit OnePointSolver( AMP::shared_ptr<SolverStrategyParameters> params ) : SolverStrategy( params )
    {
        d_onePointOp = AMP::dynamic_pointer_cast<AMP::Operator::OnePointOperator>( d_pOperator );
    }

    void solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                AMP::shared_ptr<AMP::LinearAlgebra::Vector>
                    u )
    {
        // Assumption: primaryInputVar = outputVar
        // General solution: To avoid making the above assumption, we can replace
        // getInputVariable() with getPrimaryInputVariable() and use it for the u vector.
        AMP::LinearAlgebra::Vector::shared_ptr myU =
            u->subsetVectorForVariable( d_onePointOp->getOutputVariable() );
        if ( d_bUseZeroInitialGuess ) {
            myU->zero();
        }
        AMP::LinearAlgebra::Vector::shared_ptr r = myU->cloneVector();
        d_onePointOp->residual( f, u, r );
        double inverseConstant = 1.0 / ( d_onePointOp->getConstant() );
        myU->axpy( inverseConstant, r, myU );
        // If you want to use an external solver library like Petsc or Trilinos here
        // then you will need to use the two functions defined in the
        // vectors/newFrozenVectorDesign/newFrozenVectorDesignHelpers.h file
    }

protected:
    AMP::shared_ptr<AMP::Operator::OnePointOperator> d_onePointOp;
};
}
}

#endif
