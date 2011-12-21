
#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#include "vectors/SimpleVector.h"

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

  boost::shared_ptr<AMP::Database> secondOp_db = input_db->getDatabase("FirstOperator");

  boost::shared_ptr<AMP::Database> firstOp_db = input_db->getDatabase("SecondOperator");

  boost::shared_ptr<AMP::Operator::ColumnOperator> columnOp(new AMP::Operator::ColumnOperator());

  boost::shared_ptr<AMP::Database> solver_db = input_db->getDatabase("Solver");

  boost::shared_ptr<AMP::Solver::ColumnSolver> columnSolver(new AMP::Solver::ColumnSolver());

  ut.passes(exeName);

  ut.report();
  int num_failed = ut.NumFailGlobal();

  AMP::AMPManager::shutdown();

  return num_failed;
}



