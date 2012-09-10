
#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>

#include "mpi.h"

#include "solvers/PelletStackHelpers.h"

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "ampmesh/SiloIO.h"

#include "solvers/PetscSNESSolver.h"

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

#ifdef USE_EXT_SILO
  // Create the silo writer and register the data
  AMP::Mesh::SiloIO::shared_ptr siloWriter( new AMP::Mesh::SiloIO);
#endif

  boost::shared_ptr<AMP::InputDatabase> global_input_db(new AMP::InputDatabase("global_input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, global_input_db);
  global_input_db->printClassData(AMP::plog);

  unsigned int NumberOfLoadingSteps = global_input_db->getInteger("NumberOfLoadingSteps");
  bool usePointLoad = global_input_db->getBool("USE_POINT_LOAD");
  bool useThermalLoad = global_input_db->getBool("USE_THERMAL_LOAD");

  AMP_INSIST(global_input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database> mesh_db = global_input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(globalComm);
  AMP::Mesh::Mesh::shared_ptr manager = AMP::Mesh::Mesh::buildMesh(meshParams);

  boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp;
  boost::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator;
  boost::shared_ptr<AMP::Operator::PelletStackOperator> pelletStackOp; 
  helperCreateAllOperatorsForPelletMechanics(manager, globalComm, global_input_db, coupledOp, linearColumnOperator, pelletStackOp);

  AMP::LinearAlgebra::Vector::shared_ptr solVec, rhsVec, scaledRhsVec;
  helperCreateVectorsForPelletMechanics(manager, coupledOp, solVec, rhsVec, scaledRhsVec);

  if(usePointLoad) {
    helperBuildPointLoadRHSForPelletMechanics(global_input_db, coupledOp, rhsVec);
  } else {
    rhsVec->zero();
  }

  AMP::LinearAlgebra::Vector::shared_ptr initialTemperatureVec, finalTemperatureVec;
  if(useThermalLoad) {
    helperCreateTemperatureVectorsForPelletMechanics(manager, initialTemperatureVec, finalTemperatureVec);
  }

  if(useThermalLoad) {
    double initialTemp = global_input_db->getDouble("InitialTemperature");
    initialTemperatureVec->setToScalar(initialTemp);
    helperSetReferenceTemperatureForPelletMechanics(coupledOp, initialTemperatureVec);
  }

  solVec->zero();
  helperApplyBoundaryCorrectionsForPelletMechanics(coupledOp, solVec, rhsVec);

  boost::shared_ptr<AMP::Database> nonlinearSolver_db = global_input_db->getDatabase("NonlinearSolver");
  boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver");
  boost::shared_ptr<AMP::Database> pelletStackSolver_db = linearSolver_db->getDatabase("PelletStackSolver");

  boost::shared_ptr<AMP::Solver::SolverStrategy> pelletStackSolver;
  helperBuildStackSolverForPelletMechanics(pelletStackSolver_db, pelletStackOp, linearColumnOperator, pelletStackSolver);

  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new
      AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));
  nonlinearSolverParams->d_comm = globalComm;
  nonlinearSolverParams->d_pOperator = coupledOp;
  nonlinearSolverParams->d_pInitialGuess = solVec;
  boost::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();
  linearSolver->setPreconditioner(pelletStackSolver);

#ifdef USE_EXT_SILO
  siloWriter->registerVector(solVec, manager, AMP::Mesh::Vertex, "Displacement" );
#endif

  for (unsigned int step = 0; step < NumberOfLoadingSteps; step++) {
    AMP::pout << "########################################" << std::endl;
    AMP::pout << "The current loading step is " << (step+1) << std::endl;

    double scaleValue  = (static_cast<double>(step + 1))/(static_cast<double>(NumberOfLoadingSteps));
    scaledRhsVec->scale(scaleValue, rhsVec);

    if(useThermalLoad) {
      double initialTemp = global_input_db->getDouble("InitialTemperature");
      double finalTemp = global_input_db->getDouble("FinalTemperature");
      double deltaTemp = initialTemp + ((static_cast<double>(step + 1))*(finalTemp - initialTemp)/(static_cast<double>(NumberOfLoadingSteps)));
      finalTemperatureVec->setToScalar(deltaTemp);
      helperSetFinalTemperatureForPelletMechanics(coupledOp, finalTemperatureVec);
    }

    AMP::LinearAlgebra::Vector::shared_ptr resVec = solVec->cloneVector();
    resVec->zero();
    coupledOp->apply(scaledRhsVec, solVec, resVec);
    AMP::pout<< "initial, rhsVec: "<<scaledRhsVec->L2Norm()<<endl; 
    AMP::pout<< "initial, solVec: "<<solVec->L2Norm()<<endl ; 
    AMP::pout<< "initial, resVec: "<<resVec->L2Norm()<<endl ; 
    nonlinearSolver->solve(scaledRhsVec, solVec);
    AMP::pout<< "solved,  rhsVec: "<<scaledRhsVec->L2Norm()<<endl ; 
    AMP::pout<< "solved,  solVec: "<<solVec->L2Norm()<<endl ; 
    coupledOp->apply(scaledRhsVec, solVec, resVec);
    AMP::pout<< "final,   rhsVec: "<<scaledRhsVec->L2Norm()<<endl ; 
    AMP::pout<< "final,   solVec: "<<solVec->L2Norm()<<endl ; 
    AMP::pout<< "final,   resVec: "<<resVec->L2Norm()<<endl ; 

#ifdef USE_EXT_SILO
    siloWriter->writeFile(exeName, step);
#endif

    helperResetNonlinearOperatorForPelletMechanics(coupledOp);
  }//end for step

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  int inp = 1;
  if(argc > 1) {
    inp = atoi(argv[1]);
  }

  char exeName[200];
  sprintf(exeName, "testPelletStackMechanicsSolver-%d", inp);

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


