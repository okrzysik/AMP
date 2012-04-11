
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/AMP_MPI.h"
#include "utils/PIO.h"

#include <iostream>
#include <string>

#include "vectors/MultiVariable.h"

void myTest(AMP::UnitTest *ut) {
  AMP::LinearAlgebra::Variable::shared_ptr var1(new AMP::LinearAlgebra::Variable("Hello"));
  AMP::LinearAlgebra::Variable::shared_ptr var2(new AMP::LinearAlgebra::Variable("Hello"));
  boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> multiVar( new AMP::LinearAlgebra::MultiVariable("MultiVariable"));
  multiVar->add(var1);
  multiVar->add(var2);

  AMP_ASSERT((multiVar->numVariables()) == 2);

  multiVar->removeDuplicateVariables();

  AMP_ASSERT((multiVar->numVariables()) == 1);

  ut->passes("RemDup");
}

int main(int argc, char *argv[]) {
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  try {
    myTest(&ut);
  } catch (std::exception &err) {
    AMP::pout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
    ut.failure("ERROR: While testing");
  } catch( ... ) {
    AMP::pout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
    ut.failure("ERROR: While testing");
  }

  ut.report();
  int num_failed = ut.NumFailGlobal();

  AMP::AMPManager::shutdown();
  return num_failed;
}


