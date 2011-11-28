
#include "solvers/PelletStackMechanicsSolver.h"
#include "utils/Utilities.h"

namespace AMP {
  namespace Solver {

    PelletStackMechanicsSolver :: PelletStackMechanicsSolver(boost::shared_ptr<
        SolverStrategyParameters> params) : SolverStrategy(params) {
      d_useSerial = (params->d_db)->getBoolWithDefault("USE_SERIAL", false);
      d_pelletStackOp = boost::dynamic_pointer_cast<AMP::Operator::PelletStackOperator>(d_pOperator);
    }

    void PelletStackMechanicsSolver :: setColumnSolver(boost::shared_ptr<ColumnSolver> colSolver) {
      d_columnSolver = colSolver;
    }

    void PelletStackMechanicsSolver :: resetOperator(const boost::shared_ptr<AMP::Operator::OperatorParameters> params) {
      d_columnSolver->resetOperator(params);
    }

    void PelletStackMechanicsSolver :: solve(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, 
        boost::shared_ptr<AMP::LinearAlgebra::Vector> u) {
      if(d_useSerial) {
        solveSerial(f, u);
      } else {
        solveScan(f, u);
      }
    }

    void PelletStackMechanicsSolver :: solveSerial(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, 
        boost::shared_ptr<AMP::LinearAlgebra::Vector> u) {
      AMP::LinearAlgebra::Vector::shared_ptr fCopy = f->cloneVector();

      unsigned int totalNumberOfPellets = d_pelletStackOp->getTotalNumberOfPellets();

      unsigned int locPellIdx = 0;

      if(d_pelletStackOp->hasPellet(1)) {
        boost::shared_ptr<AMP::Solver::SolverStrategy> currSolver = d_columnSolver->getSolver(locPellIdx);
        boost::shared_ptr<AMP::Operator::Operator> currOp = currSolver->getOperator();
        AMP::LinearAlgebra::Variable::shared_ptr inputVar = currOp->getInputVariable();
        AMP::LinearAlgebra::Variable::shared_ptr outputVar = currOp->getOutputVariable();
        AMP::LinearAlgebra::Vector::shared_ptr subUvec = u->subsetVectorForVariable(inputVar);
        AMP::LinearAlgebra::Vector::shared_ptr subFvec = f->subsetVectorForVariable(outputVar);
        currSolver->solve(subFvec, subUvec);
        locPellIdx++;
      }

      for(unsigned int pellId = 2; pellId <= totalNumberOfPellets; pellId++) {
        d_pelletStackOp->setCurrentPellet(pellId);
        d_pelletStackOp->apply(f, u, fCopy);
        if(d_pelletStackOp->hasPellet(pellId)) {
          boost::shared_ptr<AMP::Solver::SolverStrategy> currSolver = d_columnSolver->getSolver(locPellIdx);
          boost::shared_ptr<AMP::Operator::Operator> currOp = currSolver->getOperator();
          AMP::LinearAlgebra::Variable::shared_ptr inputVar = currOp->getInputVariable();
          AMP::LinearAlgebra::Variable::shared_ptr outputVar = currOp->getOutputVariable();
          AMP::LinearAlgebra::Vector::shared_ptr subUvec = u->subsetVectorForVariable(inputVar);
          AMP::LinearAlgebra::Vector::shared_ptr subFvec = fCopy->subsetVectorForVariable(outputVar);
          currSolver->solve(subFvec, subUvec);
          locPellIdx++;
        }
      }//end for pellId

    }

    void PelletStackMechanicsSolver :: solveScan(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, 
        boost::shared_ptr<AMP::LinearAlgebra::Vector> u) {
      AMP::LinearAlgebra::Vector::shared_ptr fCopy = f->cloneVector();
      d_columnSolver->solve(f, u);
      d_pelletStackOp->apply(f, u, fCopy);
      d_columnSolver->solve(fCopy, u);
    }

  }
}

