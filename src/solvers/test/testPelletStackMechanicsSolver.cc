
#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>

#include "mpi.h"

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "vectors/CommCollectVector.h"
#include "boost/shared_ptr.hpp"

#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshAdapter.h"
#include "ampmesh/SiloIO.h"

#include "operators/OperatorBuilder.h"
#include "operators/CoupledOperator.h"
#include "operators/map/NodeToNodeMap.h"
#include "operators/map/AsyncMapColumnOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperatorParameters.h"
#include "operators/PelletStackOperator.h"
#include "operators/boundary/DirichletVectorCorrection.h"

#include "solvers/PetscSNESSolver.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/TrilinosMLSolver.h"
#include "solvers/ColumnSolver.h"
#include "solvers/PelletStackMechanicsSolver.h"

#include "PelletStackHelpers.h"

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  boost::shared_ptr<AMP::InputDatabase> global_input_db(new AMP::InputDatabase("global_input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, global_input_db);
  global_input_db->printClassData(AMP::plog);

  unsigned int NumberOfLoadingSteps = global_input_db->getInteger("NumberOfLoadingSteps");
  bool usePointLoad = global_input_db->getBool("USE_POINT_LOAD");
  bool useThermalLoad = global_input_db->getBool("USE_THERMAL_LOAD");

  AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams ( new AMP::Mesh::MeshManagerParameters ( global_input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( mgrParams ) );

  boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps =
    AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::NodeToNodeMap> ( manager , global_input_db );

  boost::shared_ptr<AMP::Operator::PelletStackOperator> pelletStackOp; 
  helperCreatePelletStackOperator(manager, n2nmaps, globalComm, global_input_db, pelletStackOp);

  std::vector<unsigned int> localPelletIds = pelletStackOp->getLocalPelletIds();

  std::vector<AMP::Mesh::MeshManager::Adapter::shared_ptr> localMeshes = pelletStackOp->getLocalMeshes();

  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator;
  boost::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator;
  helperCreateColumnOperators(localPelletIds, localMeshes, global_input_db, nonlinearColumnOperator, linearColumnOperator);

  boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp;
  helperCreateCoupledOperator(n2nmaps, nonlinearColumnOperator, coupledOp);

  helperSetFrozenVectorForMaps(manager, globalComm, coupledOp);

  AMP::LinearAlgebra::Vector::shared_ptr solVec, rhsVec, scaledRhsVec;
  helperCreateVectors(manager, nonlinearColumnOperator, globalComm, solVec, rhsVec, scaledRhsVec);

  if(usePointLoad) {
    helperBuildPointLoadRHS(global_input_db, nonlinearColumnOperator, rhsVec);
  } else {
    rhsVec->zero();
  }

  solVec->zero();
  helperApplyBoundaryCorrections(nonlinearColumnOperator, solVec, rhsVec);

  boost::shared_ptr<AMP::Database> nonlinearSolver_db = global_input_db->getDatabase("NonlinearSolver");
  boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver");
  boost::shared_ptr<AMP::Database> pelletStackSolver_db = linearSolver_db->getDatabase("PelletStackSolver");
  boost::shared_ptr<AMP::Database> columnSolver_db = pelletStackSolver_db->getDatabase("ColumnSolver");
  boost::shared_ptr<AMP::Database> ikspSolver_db = columnSolver_db->getDatabase("KrylovSolver");
  boost::shared_ptr<AMP::Database> mlSolver_db = ikspSolver_db->getDatabase("MLSolver");

  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> columnSolverParams(new
      AMP::Solver::SolverStrategyParameters(columnSolver_db));
  columnSolverParams->d_pOperator = linearColumnOperator;
  boost::shared_ptr<AMP::Solver::ColumnSolver> columnSolver(new AMP::Solver::ColumnSolver(columnSolverParams));

  for(unsigned int id = 0; id < localPelletIds.size(); id++) {
    boost::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams(new
        AMP::Solver::SolverStrategyParameters(mlSolver_db));
    mlSolverParams->d_pOperator = linearColumnOperator->getOperator(id);
    boost::shared_ptr<AMP::Solver::TrilinosMLSolver> mlSolver(new
        AMP::Solver::TrilinosMLSolver(mlSolverParams));

    boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> ikspSolverParams(new
        AMP::Solver::PetscKrylovSolverParameters(ikspSolver_db));
    ikspSolverParams->d_pOperator = linearColumnOperator->getOperator(id);
    ikspSolverParams->d_comm = (localMeshes[id])->getComm();
    ikspSolverParams->d_pPreconditioner = mlSolver;
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> ikspSolver(new
        AMP::Solver::PetscKrylovSolver(ikspSolverParams));

    columnSolver->append(ikspSolver);
  }//end for id

  boost::shared_ptr<AMP::Solver::PelletStackMechanicsSolverParameters> pelletStackSolverParams(new
      AMP::Solver::PelletStackMechanicsSolverParameters(pelletStackSolver_db));
  pelletStackSolverParams->d_columnSolver = columnSolver;
  pelletStackSolverParams->d_pOperator = pelletStackOp;
  boost::shared_ptr<AMP::Solver::PelletStackMechanicsSolver> pelletStackSolver(new 
      AMP::Solver::PelletStackMechanicsSolver(pelletStackSolverParams));

  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new
      AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));
  nonlinearSolverParams->d_comm = globalComm;
  nonlinearSolverParams->d_pOperator = coupledOp;
  nonlinearSolverParams->d_pInitialGuess = solVec;
  boost::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));
  nonlinearSolver->setZeroInitialGuess(false);

  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();
  linearSolver->setPreconditioner(pelletStackSolver);

#ifdef USE_SILO
  manager->registerVectorAsData ( solVec , "Displacements" );
#endif

  for (unsigned int step = 0; step < NumberOfLoadingSteps; step++) {
    AMP::pout << "########################################" << std::endl;
    AMP::pout << "The current loading step is " << (step+1) << std::endl;

    double scaleValue  = ((static_cast<double>(step)) + 1.0)/(static_cast<double>(NumberOfLoadingSteps));
    scaledRhsVec->scale(scaleValue, rhsVec);

    nonlinearSolver->solve(scaledRhsVec, solVec);

#ifdef USE_SILO
    manager->writeFile<AMP::Mesh::SiloIO> ( exeName , step );
#endif

    boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
    boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperatorParameters> tmpParams(new
        AMP::Operator::MechanicsNonlinearFEOperatorParameters(tmp_db));

    for(unsigned int id = 0; id < localPelletIds.size(); id++) {
      boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinOperator =
        boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            nonlinearColumnOperator->getOperator(id));
      (nonlinOperator->getVolumeOperator())->reset(tmpParams);
    }//end for id
  }//end for step

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


