#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"
#include "materials/Material.h"


#include "ampmesh/MeshVariable.h"
#include "ampmesh/SiloIO.h"


#include "operators/MechanicsLinearFEOperator.h"
#include "operators/MechanicsNonlinearFEOperator.h"
#include "operators/DirichletVectorCorrection.h"
#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "../PetscKrylovSolverParameters.h"
#include "../PetscKrylovSolver.h"
#include "../PetscSNESSolverParameters.h"
#include "../PetscSNESSolver.h"

#include "../TrilinosMLSolver.h"


void myTest(AMP::UnitTest *ut, std::string exeName)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::AMP_MPI::initialize();  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ( "cylinder" );

  AMP_INSIST(input_db->keyExists("NumberOfLoadingSteps"), "Key ''NumberOfLoadingSteps'' is missing!");
  int NumberOfLoadingSteps = input_db->getInteger("NumberOfLoadingSteps");

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear mechanics
  AMP_INSIST( input_db->keyExists("testNonlinearMechanicsOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> mechanicsMaterialModel;
  boost::shared_ptr<AMP::Database> nonlinearMechanicsDatabase = input_db->getDatabase("testNonlinearMechanicsOperator");
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearMechanicsOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
														    "testNonlinearMechanicsOperator",
														    input_db,
														    mechanicsMaterialModel));


  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear thermal diffusion
  AMP_INSIST( input_db->keyExists("testNonlinearThermalOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
  boost::shared_ptr<AMP::Database> nonlinearThermalDatabase = input_db->getDatabase("testNonlinearThermalOperator");
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
														    "testNonlinearThermalOperator",
														    input_db,
														    thermalTransportModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a column operator object for nonlinear thermomechanics
  boost::shared_ptr<AMP::Operator::OperatorParameters> params;
  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearThermoMechanicsOperator(new AMP::Operator::ColumnOperator(params));
  nonlinearThermoMechanicsOperator->append(nonlinearMechanicsOperator);
  nonlinearThermoMechanicsOperator->append(nonlinearThermalOperator);

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the input multi-variable
  boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> volumeOperator = boost::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(nonlinearMechanicsOperator->getVolumeOperator());
  boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> inputVariable(new AMP::LinearAlgebra::MultiVariable("inputVariable"));
  inputVariable->add(volumeOperator->getInputVariable(AMP::Operator::Mechanics::DISPLACEMENT));
  inputVariable->add(volumeOperator->getInputVariable(AMP::Operator::Mechanics::TEMPERATURE));
					
  // initialize the output multi-variable
  AMP::LinearAlgebra::Variable::shared_ptr outputVariable = nonlinearThermoMechanicsOperator->getOutputVariable();

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr solVec = meshAdapter->createVector( inputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = meshAdapter->createVector( outputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = meshAdapter->createVector( outputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr workVec = meshAdapter->createVector( inputVariable );

  manager->registerVectorAsData ( solVec );
  manager->registerVectorAsData ( resVec );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // IMPORTANT:: call init before proceeding any further on the nonlinear mechanics operator
  AMP::LinearAlgebra::Vector::shared_ptr referenceTemperatureVec = meshAdapter->createVector( volumeOperator->getInputVariable(AMP::Operator::Mechanics::TEMPERATURE) );
  referenceTemperatureVec->setToScalar(300.0);
  volumeOperator->setReferenceTemperature(referenceTemperatureVec);
  volumeOperator->init();
  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // now construct the linear BVP operator for mechanics
  AMP_INSIST( input_db->keyExists("testLinearMechanicsOperator"), "key missing!" );
  boost::shared_ptr<AMP::Database> linearMechanicsDatabase = input_db->getDatabase("testLinearMechanicsOperator");
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearMechanicsOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
														 "testLinearMechanicsOperator",
														 input_db,
														 mechanicsMaterialModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // now construct the linear BVP operator for thermal
  AMP_INSIST( input_db->keyExists("testLinearThermalOperator"), "key missing!" );
  boost::shared_ptr<AMP::Database> linearThermalDatabase = input_db->getDatabase("testLinearThermalOperator");
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
														 "testLinearThermalOperator",
														 input_db,
														 thermalTransportModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a column operator object for linear thermomechanics
  boost::shared_ptr<AMP::Operator::ColumnOperator> linearThermoMechanicsOperator(new AMP::Operator::ColumnOperator(params));
  linearThermoMechanicsOperator->append(linearMechanicsOperator);
  linearThermoMechanicsOperator->append(linearThermalOperator);

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  //For RHS (Point Forces)
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Database> load_db = input_db->getDatabase("Load_Boundary");
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletLoadVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
															 "Load_Boundary",
															 input_db,
															 dummyModel));
  dirichletLoadVecOp->setVariable(residualVariable);

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  //For Initial-Guess
  boost::shared_ptr<AMP::Database> disp_db = input_db->getDatabase("Displacement_Boundary");
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletDispInVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
															 "Displacement_Boundary"
															 input_db,
															 dummyModel));
  dirichletDispInVecOp->setVariable(displacementVariable);

  
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  AMP::LinearAlgebra::Vector::shared_ptr mechNlSolVec = meshAdapter->createVector( displacementVariable );
  AMP::LinearAlgebra::Vector::shared_ptr mechNlRhsVec = meshAdapter->createVector( residualVariable );
  AMP::LinearAlgebra::Vector::shared_ptr mechNlResVec = meshAdapter->createVector( residualVariable );
  AMP::LinearAlgebra::Vector::shared_ptr mechNlScaledRhsVec = meshAdapter->createVector( residualVariable );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  //Initial guess for NL solver must satisfy the displacement boundary conditions
  mechNlSolVec->setRandomValues();
  dirichletDispInVecOp->apply(nullVec, nullVec, mechNlSolVec, 1.0, 0.0);

  nonlinearMechanicsOperator->apply(nullVec, mechNlSolVec, mechNlResVec, 1.0, 0.0);
  linearMechanicsOperator->reset(nonlinBvpOperator->getJacobianParameters(mechNlSolVec));

  mechNlRhsVec->setToScalar(0.0);
  dirichletLoadVecOp->apply(nullVec, nullVec, mechNlRhsVec, 1.0, 0.0);

  double initSolNorm = mechNlSolVec->L2Norm();

  std::cout<<"Initial Solution Norm: "<<initSolNorm<<std::endl;

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 

  boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // ---- first initialize the preconditioner for mechanics
  boost::shared_ptr<AMP::Database> mechanicsPCSolver_db = linearSolver_db->getDatabase("MechanicsPreconditioner"); 
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> mechanicsPCSolverParams(new AMP::Solver::SolverStrategyParameters(mechanicsPCSolver_db));
  mechanicsPCSolverParams->d_pOperator = linearMechanicsOperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> mechanicsPCSolver(new AMP::Solver::TrilinosMLSolver(pcSolverParams));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // ----  initialize the preconditioner for thermal
  boost::shared_ptr<AMP::Database> thermalPCSolver_db = linearSolver_db->getDatabase("ThermalPreconditioner"); 
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPCSolverParams(new AMP::Solver::SolverStrategyParameters(thermalPCSolver_db));
  thermalPCSolverParams->d_pOperator = linearThermalOperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> thermalPCSolver(new AMP::Solver::TrilinosMLSolver(thermalPCSolverParams));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // construct the composite preconditioner
  
  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the linear solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));

  linearSolverParams->d_pOperator = linearThermoMechanicsOperator;
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pPreconditioner = pcSolver;

  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the nonlinear solver
  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new
      AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm = globalComm;
  nonlinearSolverParams->d_pOperator = nonlinearThermoMechanicsOperator;
  nonlinearSolverParams->d_pKrylovSolver = linearSolver;
  nonlinearSolverParams->d_pInitialGuess = nonlinearSolVec;

  boost::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

  nonlinearSolver->setZeroInitialGuess(false);

  linearThermoMechanicsOperator->reset(nonlinearThermoMechanicsOperator->getJacobianParameters(nonlinearSolVec));

  for (int step=0;step<NumberOfLoadingSteps; step++)
  {
    std::cout << "########################################" << std::endl;
    std::cout << "The current loading step is " << (step+1) << std::endl;
    double scaleValue  = ((double)step+1.0)/NumberOfLoadingSteps;
    mechNlScaledRhsVec->scale(scaleValue, mechNlRhsVec);
    std:: cout << "L2 Norm at loading step " << (step+1) << " is " << mechNlScaledRhsVec->L2Norm() << std::endl;
    nonlinearSolver->solve(mechNlScaledRhsVec, mechNlSolVec);
    boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
    boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperatorParameters> tmpParams(new
        AMP::Operator::MechanicsNonlinearFEOperatorParameters(tmp_db));
    nonlinearThermoMechanicsOperator->reset(tmpParams);
    nonlinearSolver->setZeroInitialGuess(false);
  }

  double finalSolNorm = mechNlSolVec->L2Norm();

  std::cout<<"Final Solution Norm: "<<finalSolNorm<<std::endl;

  nonlinBvpOperator->apply(mechNlRhsVec, mechNlSolVec, mechNlResVec, 1.0, -1.0);

  double finalResidualNorm  = mechNlResVec->L2Norm();

  std::cout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;

#ifdef USE_SILO
  manager->writeFile<AMP::Mesh::SiloIO> ( "vonMisesExample" , 1 );
#endif

  if(finalResidualNorm>initialResidualNorm*1.0e-10+1.0e-05)
  {
    ITFAILS;
  }
  else
  {
    ut.passes("PetscSNES Solver successfully solves a nonlinear mechanics equation with Jacobian provided, FGMRES for Krylov");
  }

  ut.passes(exeName);

  AMP::AMPManager::shutdown();

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

  node  = AMP::AMP_MPI::getRank();
  nodes = AMP::AMP_MPI::getNodes();

  std::vector<std::string> exeNames;
  exeNames.push_back("testPetscSNESSolver-ML-NonlinearMechanics-1-normal");
  exeNames.push_back("testPetscSNESSolver-ML-NonlinearMechanics-1-reduced");

  for(int i = 0; i < exeNames.size(); i++) {
    try {      myTest(ut, exeNames[i]);    }
    catch (std::exception &err)
    {
      std::cout << "ERROR: While testing "<<argv[0] 
        << err.what()
        << std::endl;
    }
    catch( ... )
    {
      std::cout << "ERROR: While testing "<<argv[0] 
        << "An unknown exception was thrown."
        << std::endl;
    }
  } //end for i

  return ut.numFails;
}   


