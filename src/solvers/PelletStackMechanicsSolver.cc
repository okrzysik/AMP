
#include "solvers/PelletStackMechanicsSolver.h"

namespace AMP {
  namespace Solver {

    PelletStackMechanicsSolver :: PelletStackMechanicsSolver(boost::shared_ptr<
        PelletStackMechanicsSolverParameters> params) : SolverStrategy(params) {
      d_columnSolver = params->d_columnSolver;
      d_pelletStackOp = boost::dynamic_pointer_cast<AMP::Operator::PelletStackOperator>(d_pOperator);
      d_fCopy.reset();
    }

    void PelletStackMechanicsSolver :: resetOperator(const boost::shared_ptr<AMP::Operator::OperatorParameters> params) {
      d_columnSolver->resetOperator(params);
    }

    void PelletStackMechanicsSolver :: solve(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, 
        boost::shared_ptr<AMP::LinearAlgebra::Vector> u) {
      boost::shared_ptr<AMP::LinearAlgebra::Vector> fInternal = f;
      if(d_pelletStackOp->useScaling()) {
        if(d_fCopy == NULL) {
          d_fCopy = f->cloneVector();
        }
        d_fCopy->copyVector(f);
        d_pelletStackOp->applyUnscaling(d_fCopy);
        fInternal = d_fCopy;
      }
      if(d_pelletStackOp->useSerial()) {
        solveSerial(fInternal, u);
      } else {
        solveScan(fInternal, u);
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
      d_columnSolver->solve(f, u);
      if(d_pelletStackOp->onlyZcorrection()) {
        AMP::LinearAlgebra::Vector::shared_ptr nullVec;
        d_pelletStackOp->apply(nullVec, nullVec, u);
      } else {
        AMP::LinearAlgebra::Vector::shared_ptr fCopy = f->cloneVector();
        d_pelletStackOp->apply(f, u, fCopy);
        d_columnSolver->solve(fCopy, u);
      }
    }

  }
}


