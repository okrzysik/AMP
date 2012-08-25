#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/PIO.h"
#include "materials/Material.h"


#include "ampmesh/MeshVariable.h"
#include "ampmesh/SiloIO.h"


#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/ColumnOperator.h"
#include "operators/OperatorBuilder.h"

#include "../ColumnSolver.h"
#include "../PetscKrylovSolverParameters.h"
#include "../PetscKrylovSolver.h"
#include "../NonlinearKrylovAcceleratorParameters.h"
#include "../NonlinearKrylovAccelerator.h"

#include "../TrilinosMLSolver.h"


void myTest(AMP::UnitTest *ut, std::string exeName)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST(input_db->keyExists("NumberOfMeshes"), "Key does not exist");
  int numMeshes = input_db->getInteger("NumberOfMeshes");
  
  AMP::pout<<"Num meshes = "<<numMeshes<<std::endl;

  AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ( "cylinder" );
  
  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear thermal diffusion
  AMP_INSIST( input_db->keyExists("testNonlinearThermalOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
  boost::shared_ptr<AMP::Database> nonlinearThermalDatabase = input_db->getDatabase("testNonlinearThermalOperator");
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator = boost::dynamic_pointer_cast<
    AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
											"testNonlinearThermalOperator",
											input_db,
											thermalTransportModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the input variable
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator =
		  boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nonlinearThermalOperator->getVolumeOperator());

  boost::shared_ptr<AMP::LinearAlgebra::Variable> inputVariable(thermalVolumeOperator->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));

  // initialize the output variable
  AMP::LinearAlgebra::Variable::shared_ptr outputVariable = thermalVolumeOperator->getOutputVariable();

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr solVec = meshAdapter->createVector( inputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = meshAdapter->createVector( outputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = meshAdapter->createVector( outputVariable );

  // create the following shared pointers for ease of use
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // register some variables for plotting
  meshAdapter->registerVectorAsData ( solVec, "Solution" );
  meshAdapter->registerVectorAsData ( resVec, "Residual" );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // now construct the linear BVP operator for thermal
  AMP_INSIST( input_db->keyExists("testLinearThermalOperator"), "key missing!" );
  boost::shared_ptr<AMP::Database> linearThermalDatabase = input_db->getDatabase("testLinearThermalOperator");
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
																	    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																							   "testLinearThermalOperator",
																							   input_db,
																							   thermalTransportModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  //Initial-Guess for thermal
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyThermalModel;
  boost::shared_ptr<AMP::Database> thermal_db = input_db->getDatabase("ThermalInitialGuess");
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletThermalInVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
															 "ThermalInitialGuess",
															 input_db,
															 dummyThermalModel));
  dirichletThermalInVecOp->setVariable(thermalVolumeOperator->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  //Initial guess

  solVec->setToScalar(400.);
  double initialGuessNorm  = solVec->L2Norm();
  std::cout << "initial guess norm = " << initialGuessNorm <<"\n";

  //Initial guess for thermal must satisfy the thermal Dirichlet boundary conditions
  dirichletThermalInVecOp->apply(nullVec, nullVec, solVec, 1.0, 0.0);

  initialGuessNorm  = solVec->L2Norm();
  std::cout << "initial guess norm  after apply = " << initialGuessNorm <<"\n";

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  //Enforce consistent rhs
  boost::shared_ptr<AMP::Operator::BoundaryOperator> thermalBoundaryOperator =
		  nonlinearThermalOperator->getBoundaryOperator();

  rhsVec->setToScalar(0.5);

  thermalBoundaryOperator->apply(nullVec, nullVec, rhsVec, 1.0, 0.0);

  double rhsNorm  = rhsVec->L2Norm();
  std::cout << "rhs norm  after apply = " << rhsNorm <<"\n";

  linearThermalOperator->reset(nonlinearThermalOperator->getJacobianParameters(solVec));
  
  //----------------------------------------------------------------------------------------------------------------------------------------------/

  boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
  boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  boost::shared_ptr<AMP::Database> thermalPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams(new AMP::Solver::SolverStrategyParameters(thermalPreconditioner_db));
  thermalPreconditionerParams->d_pOperator = linearThermalOperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> linearThermalPreconditioner(new AMP::Solver::TrilinosMLSolver(thermalPreconditionerParams));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the nonlinear solver
  boost::shared_ptr<AMP::Solver::NonlinearKrylovAcceleratorParameters> nonlinearSolverParams(new
      AMP::Solver::NonlinearKrylovAcceleratorParameters(nonlinearSolver_db));

  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_pOperator = nonlinearThermalOperator;
  nonlinearSolverParams->d_pInitialGuess = solVec;
  nonlinearSolverParams->d_pPreconditioner = linearThermalPreconditioner;
  boost::shared_ptr<AMP::Solver::NonlinearKrylovAccelerator> nonlinearSolver(new AMP::Solver::NonlinearKrylovAccelerator(nonlinearSolverParams));


  nonlinearThermalOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);
  double initialResidualNorm  = resVec->L2Norm();

  AMP::pout<<"Initial Residual Norm: "<<initialResidualNorm<<std::endl;

  nonlinearSolver->setZeroInitialGuess(false);

  nonlinearSolver->solve(rhsVec, solVec);

  nonlinearThermalOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);

  double finalResidualNorm  = resVec->L2Norm();

  std::cout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;

#ifdef USES_SILO
  manager->writeFile<AMP::Mesh::SiloIO> ( exeName , 0 );
#endif

  if(finalResidualNorm>1.0e-08) {
    ut->failure("Error");
  } else {
    ut->passes("NonlinearKrylovAccelerator successfully solves a nonlinear thermal equation with Jacobian provided, FGMRES for Krylov");
  }
  ut->passes(exeName);

}


int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back("testNonlinearKrylovAccelerator-NonlinearThermal-cylinder_kIsOne");
    exeNames.push_back("testNonlinearKrylovAccelerator-NonlinearThermal-cylinder_MATPRO");

    for(unsigned int i = 0; i < exeNames.size(); i++) {
        try {
            myTest(&ut, exeNames[i]);
        } catch (std::exception &err) {
            std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
            ut.failure("ERROR: While testing");
        } catch( ... ) {
            std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
            ut.failure("ERROR: While testing");
        }
    }
   
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   



