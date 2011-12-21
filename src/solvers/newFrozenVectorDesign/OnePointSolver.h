
#ifndef included_AMP_OnePointSolver
#define included_AMP_OnePointSolver

#include "solvers/SolverStrategy.h"
#include "operators/newFrozenVectorDesign/OnePointOperator.h"
#include "vectors/newFrozenVectorDesign/newFrozenVectorDesignHelpers.h"

namespace AMP {
  namespace Solver {

    class OnePointSolver: public SolverStrategy {
      public:
        OnePointSolver(boost::shared_ptr<SolverStrategyParameters> params) : SolverStrategy(params) {
          d_onePointOp = boost::dynamic_pointer_cast<AMP::Operator::OnePointOperator>(d_pOperator);
        }

        void solve(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, boost::shared_ptr<AMP::LinearAlgebra::Vector> u) {
          //Assumption: primaryInputVar = outputVar
          //General solution: To avoid making the above assumption, we can replace 
          //getInputVariable() with getPrimaryVariable() and use it for the u vector. 
          AMP::LinearAlgebra::Vector::shared_ptr myU = u->subsetVectorForVariable(d_onePointOp->getOutputVariable());
          AMP::LinearAlgebra::Vector::shared_ptr myF = f->subsetVectorForVariable(d_onePointOp->getOutputVariable());
          if(d_bUseZeroInitialGuess) {
            myU->zero();
          }
          AMP::LinearAlgebra::Vector::shared_ptr r = myU->cloneVector();
          d_onePointOp->apply(f, u, r, -1.0, 1.0);
          double inverseConstant = 1.0/(d_onePointOp->getConstant());
          myU->axpy(inverseConstant, r, myU);
          //If you want to use an external solver library like Petsc or Trilinos here 
          //then you will need to use the two functions defined in the
          //vectors/newFrozenVectorDesign/newFrozenVectorDesignHelpers.h file
        }

      protected:
        boost::shared_ptr<AMP::Operator::OnePointOperator> d_onePointOp;
    };

  }
}

#endif



