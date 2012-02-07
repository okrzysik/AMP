#include <string>
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "materials/Material.h"
#include "boost/shared_ptr.hpp"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "operators/NeutronicsRhs.h"
#include "vectors/Variable.h"

#include "ampmesh/SiloIO.h"
#include "vectors/Vector.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/VolumeIntegralOperator.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/ElementOperationFactory.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/boundary/RobinMatrixCorrection.h"
#include "operators/boundary/NeumannVectorCorrection.h"

#include "../TrilinosMLSolver.h"


void linearThermalTest(AMP::UnitTest *ut )
{
  //double t1;
  // Input and output file names
  //  #include <string>
  std::string exeName("testTrilinosMLSolver-LinearThermalRobin");
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;
////////////////////////////////////
//    INITIALIZE THE PROBLEM      //
////////////////////////////////////

  // Construct a smart pointer to a new database.
  //  #include "boost/shared_ptr.hpp"
  //  #include "utils/InputDatabase.h"
  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));

  // Fill the database from the input file.
  //  #include "utils/InputManager.h"
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);


  // Print from all cores into the output files
  //   #include "utils/PIO.h"
  AMP::PIO::logAllNodes(log_file);

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  //std::string mesh_file = input_db->getString("Mesh");

  // Construct a mesh manager which reads in the fuel mesh
  AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( mgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ( "bar" );
  
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
////////////////////////////////////
//  CREATE THE NEUTRONICS SOURCE  //
////////////////////////////////////

  AMP_INSIST(input_db->keyExists("NeutronicsOperator"), "Key ''NeutronicsOperator'' is missing!");
  boost::shared_ptr<AMP::Database>  neutronicsOp_db = input_db->getDatabase("NeutronicsOperator");
  boost::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ));
  neutronicsParams->d_MeshAdapter = meshAdapter;
  boost::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(new AMP::Operator::NeutronicsRhs( neutronicsParams ));

  AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar = neutronicsOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   SpecificPowerVec = meshAdapter->createVector( SpecificPowerVar );

  neutronicsOperator->apply(nullVec, nullVec, SpecificPowerVec, 1., 0.);

/////////////////////////////////////////////////////
//  Integrate Nuclear Source over Desnity * Volume //
/////////////////////////////////////////////////////

  AMP_INSIST( input_db->keyExists("VolumeIntegralOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
  boost::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator = boost::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																							      "VolumeIntegralOperator",
																							      input_db,
																							      stransportModel));

  // Create the power (heat source) vector.
  AMP::LinearAlgebra::Variable::shared_ptr PowerInWattsVar = sourceOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   PowerInWattsVec = meshAdapter->createVector( PowerInWattsVar );
  PowerInWattsVec->zero();

  // convert the vector of specific power to power for a given basis.
  sourceOperator->apply(nullVec, SpecificPowerVec, PowerInWattsVec, 1., 0.);

////////////////////////////////////////
//   CREATE THE THERMAL BVP OPERATOR  //
////////////////////////////////////////
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> diffusionOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																						       "DiffusionBVPOperator",
																						       input_db,
																						       transportModel));


  AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec = meshAdapter->createVector( diffusionOperator->getInputVariable() );
  AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec       = meshAdapter->createVector( diffusionOperator->getOutputVariable() );
  AMP::LinearAlgebra::Vector::shared_ptr ResidualVec            = meshAdapter->createVector( diffusionOperator->getOutputVariable() );

  RightHandSideVec->setToScalar(0.0);
  double rhsNorm = RightHandSideVec->L2Norm();

///////////////////////////////////////////////
//   Add the boundary conditions corrections //
///////////////////////////////////////////////

  AMP::LinearAlgebra::Vector::shared_ptr boundaryOpCorrectionVec = meshAdapter->createVector( diffusionOperator->getOutputVariable() );

  AMP::Operator::Operator::shared_ptr  boundaryOp;
  boundaryOp = diffusionOperator->getBoundaryOperator();

  (boost::dynamic_pointer_cast<AMP::Operator::BoundaryOperator>(boundaryOp) )->addRHScorrection(boundaryOpCorrectionVec) ;

  RightHandSideVec->subtract(PowerInWattsVec, boundaryOpCorrectionVec);

  rhsNorm = RightHandSideVec->L2Norm();
  std::cout<<"RHS Norm after BC Correction "<<rhsNorm<<std::endl;
 
  rhsNorm = RightHandSideVec->L2Norm();
  std::cout<<"RHS Norm 1: "<<rhsNorm<<std::endl;
  rhsNorm = PowerInWattsVec->L2Norm();
  std::cout<<"RHS Norm 2: "<<rhsNorm<<std::endl;
  rhsNorm = boundaryOpCorrectionVec->L2Norm();
  std::cout<<"RHS Norm 3: "<<rhsNorm<<std::endl;

/////////////////////////////////////////////
  // make sure the database on theinput file exists for the linear solver
  AMP_INSIST(input_db->keyExists("LinearSolver"),   "Key ''LinearSolver'' is missing!");

  // Read the input file onto a database.
  boost::shared_ptr<AMP::Database>                 mlSolver_db   = input_db->getDatabase("LinearSolver"); 

  // Fill in the parameters fo the class with the info on the database.
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams (new AMP::Solver::SolverStrategyParameters(mlSolver_db));

  // Define the operature to be used by the Solver.
  mlSolverParams->d_pOperator = diffusionOperator;

  // Set initial guess
  TemperatureInKelvinVec->setToScalar(1.0);

  // Check the initial L2 norm of the solution
  double initSolNorm = TemperatureInKelvinVec->L2Norm();
  std::cout<<"Initial Solution Norm: "<<initSolNorm<<std::endl;

  rhsNorm = RightHandSideVec->L2Norm();
  std::cout<<"RHS Norm: "<<rhsNorm<<std::endl;

  // Create the ML Solver
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver>         mlSolver(new AMP::Solver::TrilinosMLSolver(mlSolverParams));
 
  // Use a random initial guess?
  mlSolver->setZeroInitialGuess(false);

  // Solve the prblem.
  mlSolver->solve(RightHandSideVec, TemperatureInKelvinVec);

  // Compute the residual
  diffusionOperator->apply(RightHandSideVec, TemperatureInKelvinVec, ResidualVec);

  // Check the L2 norm of the final residual.
  double finalResidualNorm = ResidualVec->L2Norm();
  std::cout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;

  if(finalResidualNorm>10.0) {
    ut->failure("TrilinosMLSolver successfully solves a linear thermal problem with a nuclear source term.");
  } else {
    ut->passes("TrilinosMLSolver successfully solves a linear thermal problem with a nuclear source term.");
  }

  // Plot the results
  AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
  if( globalComm.getSize() == 1 ) {
#ifdef USE_SILO
      AMP::LinearAlgebra::Variable::shared_ptr tmpVar1 = TemperatureInKelvinVec->getVariable();
      tmpVar1->setName("TemperatureInKelvin");
      
      tmpVar1 = ResidualVec->getVariable();
      tmpVar1->setName("Residual");
      
//    manager->getMesh ( "fuel" )->registerVectorAsData ( SpecificPowerVec );
/////////    manager->getMesh ( "fuel" )->registerVectorAsData ( ResidualVec );
      meshAdapter->registerVectorAsData ( ResidualVec );
//    manager->getMesh ( "fuel" )->registerVectorAsData ( PowerInWattsVec );
/////////    manager->getMesh ( "fuel" )->registerVectorAsData ( TemperatureInKelvinVec );
      meshAdapter->registerVectorAsData ( TemperatureInKelvinVec );
//    manager->getMesh ( "fuel" )->registerVectorAsData ( RightHandSideVec );
    manager->writeFile<AMP::Mesh::SiloIO> ( exeName, 0 );
#endif
  }

  input_db.reset();

  ut->passes(exeName);

}


int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
        linearThermalTest(&ut);
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



