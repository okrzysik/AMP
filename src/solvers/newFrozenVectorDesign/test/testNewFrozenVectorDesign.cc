
#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#include "vectors/SimpleVector.h"
#include "vectors/MultiVector.h"
#include "vectors/newFrozenVectorDesign/newFrozenVectorDesignHelpers.h"

#include "operators/ColumnOperator.h"
#include "operators/newFrozenVectorDesign/FirstOperator.h"
#include "operators/newFrozenVectorDesign/SecondOperator.h"

#include "solvers/newFrozenVectorDesign/OnePointSolver.h"
#include "solvers/ColumnSolver.h"

int main(int argc, char *argv[]) {
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testNewFrozenVectorDesign";
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;
  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);

  boost::shared_ptr<AMP::Database> firstOp_db = input_db->getDatabase("FirstOperator");
  boost::shared_ptr<AMP::Operator::OperatorParameters> firstOpParams(new AMP::Operator::OperatorParameters(firstOp_db));
  boost::shared_ptr<AMP::Operator::FirstOperator> firstOp(new AMP::Operator::FirstOperator(firstOpParams)); 

  boost::shared_ptr<AMP::Database> secondOp_db = input_db->getDatabase("SecondOperator");
  boost::shared_ptr<AMP::Operator::OperatorParameters> secondOpParams(new AMP::Operator::OperatorParameters(secondOp_db));
  boost::shared_ptr<AMP::Operator::SecondOperator> secondOp(new AMP::Operator::SecondOperator(secondOpParams)); 

  boost::shared_ptr<AMP::Operator::ColumnOperator> columnOp(new AMP::Operator::ColumnOperator);
  columnOp->append(firstOp);
  columnOp->append(secondOp);

  std::cout<<"Constructed Operators."<<std::endl;

  boost::shared_ptr<AMP::Database> solver_db = input_db->getDatabase("Solver");

  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> firstSolverParams(new AMP::Solver::SolverStrategyParameters(solver_db));
  firstSolverParams->d_pOperator = firstOp;
  boost::shared_ptr<AMP::Solver::OnePointSolver> firstSolver(new AMP::Solver::OnePointSolver(firstSolverParams));

  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> secondSolverParams(new AMP::Solver::SolverStrategyParameters(solver_db));
  secondSolverParams->d_pOperator = secondOp;
  boost::shared_ptr<AMP::Solver::OnePointSolver> secondSolver(new AMP::Solver::OnePointSolver(secondSolverParams));

  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> columnSolverParams(new AMP::Solver::SolverStrategyParameters(solver_db));
  columnSolverParams->d_pOperator = columnOp;
  boost::shared_ptr<AMP::Solver::ColumnSolver> columnSolver(new AMP::Solver::ColumnSolver(columnSolverParams));
  columnSolver->append(firstSolver);
  columnSolver->append(secondSolver);

  std::cout<<"Constructed Solvers."<<std::endl;

  AMP::LinearAlgebra::Variable::shared_ptr firstVar = firstOp->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr secondVar = secondOp->getOutputVariable();

  AMP::LinearAlgebra::Vector::shared_ptr firstVecTmp = AMP::LinearAlgebra::SimpleVector::create(1, firstVar);
  AMP::LinearAlgebra::Vector::shared_ptr secondVecTmp = AMP::LinearAlgebra::SimpleVector::create(1, secondVar);

  AMP::LinearAlgebra::Vector::shared_ptr firstVec = AMP::LinearAlgebra::MultiVector::view(firstVecTmp, AMP_COMM_SELF);
  AMP::LinearAlgebra::Vector::shared_ptr secondVec = AMP::LinearAlgebra::MultiVector::view(secondVecTmp, AMP_COMM_SELF);

  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::joinVectors(firstVec, secondVec);
  AMP::LinearAlgebra::Vector::shared_ptr solVec = rhsVec->cloneVector();

  std::cout<<"Constructed Vectors."<<std::endl;

  firstVec->setToScalar(3.0);
  secondVec->setToScalar(5.0);

  std::cout<<"Formed RHS."<<std::endl;

  columnSolver->solve(rhsVec, solVec);

  std::cout<<"Completed Solve."<<std::endl;

  boost::shared_ptr<AMP::LinearAlgebra::SimpleVector> firstSol = 
    boost::dynamic_pointer_cast<AMP::LinearAlgebra::SimpleVector>(solVec->subsetVectorForVariable(firstVar));

  boost::shared_ptr<AMP::LinearAlgebra::SimpleVector> secondSol = 
    boost::dynamic_pointer_cast<AMP::LinearAlgebra::SimpleVector>(solVec->subsetVectorForVariable(secondVar));

  std::cout<<"First Solution = "<<((*firstSol)[0])<<std::endl;
  std::cout<<"Second Solution = "<<((*secondSol)[0])<<std::endl;

  ut.passes(exeName);

  ut.report();
  int num_failed = ut.NumFailGlobal();

  AMP::AMPManager::shutdown();

  return num_failed;
}



