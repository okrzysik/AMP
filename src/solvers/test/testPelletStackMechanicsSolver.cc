
#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "mpi.h"

#include "vectors/CommCollectVector.h"
#include "boost/shared_ptr.hpp"

#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshAdapter.h"

#include "operators/OperatorBuilder.h"
#include "operators/CoupledOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/PelletStackOperator.h"
#include "operators/boundary/DirichletVectorCorrection.h"

#include "solvers/PetscSNESSolver.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/TrilinosMLSolver.h"
#include "solvers/ColumnSolver.h"
#include "solvers/PelletStackMechanicsSolver.h"

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  boost::shared_ptr<AMP::InputDatabase> global_input_db(new AMP::InputDatabase("global_input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, global_input_db);
  global_input_db->printClassData(AMP::plog);
  
  unsigned int NumberOfLoadingSteps = global_input_db->getInteger("NumberOfLoadingSteps");
  bool usePointForceLoad = global_input_db->getBool("USE_POINT_FORCE_LOAD");
  bool useThermalLoad = global_input_db->getBool("USE_THERMAL_LOAD");
  
  AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams ( new AMP::Mesh::MeshManagerParameters ( global_input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( mgrParams ) );

  boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  n2nmaps =
    AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::NodeToNodeMap> ( manager , global_input_db ) ;

  boost::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator(new AMP::Operator::ColumnOperator(emptyParams));
  boost::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator(new AMP::Operator::ColumnOperator(emptyParams));

  boost::shared_ptr<AMP::InputDatabase> emptyDb;
  boost::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledOpParams(new
      AMP::Operator::CoupledOperatorParameters(emptyDb));
  coupledOpParams->d_MapOperator = n2nmaps;
  coupledOpParams->d_BVPOperator = nonlinearColumnOperator;
  boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp(new AMP::Operator::CoupledOperator(coupledOpParams));

  std::vector<boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> > pointLoadOperators;

  std::vector<unsigned int> localMeshIds;

  std::vector<AMP::Mesh::MeshManager::Adapter::shared_ptr> localMeshes;

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testPelletStackMechanicsSolver";

  try {
    myTest(&ut, exeName);
  } catch (std::exception &err) {
    std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
    ut.failure("ERROR: While testing");
  } catch( ... ) {
    std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
    ut.failure("ERROR: While testing");
  }

  ut.report();

  int num_failed = ut.NumFailGlobal();
  AMP::AMPManager::shutdown();
  return num_failed;
}  


