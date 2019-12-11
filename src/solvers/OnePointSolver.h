
#ifndef included_AMP_OnePointSolver
#define included_AMP_OnePointSolver

#include "AMP/operators/newFrozenVectorDesign/OnePointOperator.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/vectors/newFrozenVectorDesign/newFrozenVectorDesignHelpers.h"

namespace AMP {
namespace Solver {

class OnePointSolver : public SolverStrategy
{
public:
    explicit OnePointSolver( std::shared_ptr<SolverStrategyParameters> params )
        : SolverStrategy( params ),
          d_onePointOp( std::dynamic_pointer_cast<AMP::Operator::OnePointOperator>( d_pOperator ) )
    {
    }

    virtual void solve( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                        std::shared_ptr<AMP::LinearAlgebra::Vector> u ) override
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
    std::shared_ptr<AMP::Operator::OnePointOperator> d_onePointOp;
};
} // namespace Solver
} // namespace AMP

#endif
