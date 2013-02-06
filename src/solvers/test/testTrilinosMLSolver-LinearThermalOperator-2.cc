
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

#include "ampmesh/Mesh.h"
#include "vectors/VectorBuilder.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "utils/Writer.h"
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
#include "operators/boundary/NeumannVectorCorrection.h"

#include "solvers/trilinos/TrilinosMLSolver.h"


void linearThermalTest(AMP::UnitTest *ut )
{
  double t1;
  // Input and output file names
  //  #include <string>
  std::string exeName("testTrilinosMLSolver-LinearThermalOperator-2_HALDEN");
  std::string input_file = "input_" + exeName;
  //std::string input_file = "input_" + exeName + "_HALDEN";
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

//--------------------------------------------------
//   Create the Mesh.
//--------------------------------------------------
    AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
    boost::shared_ptr<AMP::Database>  mesh_db = input_db->getDatabase("Mesh");
    boost::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(new AMP::Mesh::MeshParameters(mesh_db));
    mgrParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
    boost::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh(mgrParams);
//--------------------------------------------------

//--------------------------------------------------
// Create a DOF manager for a nodal vector 
//--------------------------------------------------
  int DOFsPerNode = 1;
  int DOFsPerElement = 8;
  int nodalGhostWidth = 1;
  int gaussPointGhostWidth = 1;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr nodalDofMap      = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth,      DOFsPerNode,    split);
  AMP::Discretization::DOFManager::shared_ptr gaussPointDofMap = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Volume, gaussPointGhostWidth, DOFsPerElement, split);
//--------------------------------------------------

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  ////////////////////////////////////
  //  CREATE THE NEUTRONICS SOURCE  //
  ////////////////////////////////////
  AMP_INSIST(input_db->keyExists("NeutronicsOperator"), "Key ''NeutronicsOperator'' is missing!");
  boost::shared_ptr<AMP::Database>  neutronicsOp_db = input_db->getDatabase("NeutronicsOperator");
  boost::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ));
  boost::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(new AMP::Operator::NeutronicsRhs( neutronicsParams ));

  AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar = neutronicsOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   SpecificPowerVec = AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );

  neutronicsOperator->apply(nullVec, nullVec, SpecificPowerVec, 1., 0.);

  /////////////////////////////////////////////////////
  //  Integrate Nuclear Rhs over Desnity * Volume //
  /////////////////////////////////////////////////////

  AMP_INSIST( input_db->keyExists("VolumeIntegralOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
  boost::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator = boost::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
																	       AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																							      "VolumeIntegralOperator",
																							      input_db,
																							      stransportModel));

  // Create the power (heat source) vector.
  AMP::LinearAlgebra::Variable::shared_ptr PowerInWattsVar = sourceOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   PowerInWattsVec = AMP::LinearAlgebra::createVector( nodalDofMap, PowerInWattsVar );
  PowerInWattsVec->zero();

  // convert the vector of specific power to power for a given basis.
  sourceOperator->apply(nullVec, SpecificPowerVec, PowerInWattsVec, 1., 0.);

  t1 = SpecificPowerVec->L2Norm();
  std::cout << "n1 = " << t1 << std::endl;
  t1 = PowerInWattsVec->L2Norm();
  std::cout << "n1 = " << t1 << std::endl;

  ////////////////////////////////////
  //   CREATE THE THERMAL OPERATOR  //
  ////////////////////////////////////
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> diffusionOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
																	AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																						       "DiffusionBVPOperator",
																						       input_db,
																						       transportModel));

  AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec = AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getInputVariable()  );
  AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec       = AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable() );
  AMP::LinearAlgebra::Vector::shared_ptr ResidualVec            = AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable() );

  RightHandSideVec->setToScalar(0.0);
  double rhsNorm = RightHandSideVec->L2Norm();

  ///////////////////////////////////////////////
  //   Add the boundary conditions corrections //
  ///////////////////////////////////////////////

  RightHandSideVec->copyVector(PowerInWattsVec);

  diffusionOperator->modifyRHSvector(RightHandSideVec);

  rhsNorm = RightHandSideVec->L2Norm();
  std::cout<<"RHS Norm 1: "<<rhsNorm<<std::endl;
  rhsNorm = PowerInWattsVec->L2Norm();
  std::cout<<"RHS Norm 2: "<<rhsNorm<<std::endl;


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

  AMP::pout<<"RHS Max: "<<RightHandSideVec->max()<<std::endl;
  AMP::pout<<"RHS Min: "<<RightHandSideVec->min()<<std::endl;
  AMP::pout<<"RHS L1-norm: "<<RightHandSideVec->L1Norm()<<std::endl;
  AMP::pout<<"RHS L2-norm: "<<RightHandSideVec->L2Norm()<<std::endl;
  
  // Solve the prblem.
  mlSolver->solve(RightHandSideVec, TemperatureInKelvinVec);

  AMP::pout<<"Solution Max: "<<TemperatureInKelvinVec->max()<<std::endl;
  AMP::pout<<"Solution Min: "<<TemperatureInKelvinVec->min()<<std::endl;
  AMP::pout<<"Solution L1-norm: "<<TemperatureInKelvinVec->L1Norm()<<std::endl;
  AMP::pout<<"Solution L2-norm: "<<TemperatureInKelvinVec->L2Norm()<<std::endl;
  
  // Compute the residual
  diffusionOperator->apply(RightHandSideVec, TemperatureInKelvinVec, ResidualVec);

  // Check the L2 norm of the final residual.
  double finalResidualNorm = ResidualVec->L2Norm();
  AMP::pout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;

  if(finalResidualNorm>10.0) {
    ut->failure("TrilinosMLSolver successfully solves a linear thermal problem with a nuclear source term.");
  } else {
    ut->passes("TrilinosMLSolver successfully solves a linear thermal problem with a nuclear source term.");
  }

  // Plot the results
  AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
#ifdef USE_EXT_SILO
     AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
     siloWriter->registerMesh( meshAdapter );

     siloWriter->registerVector( TemperatureInKelvinVec, meshAdapter, AMP::Mesh::Vertex, "TemperatureInKelvin" );
     siloWriter->registerVector( ResidualVec,            meshAdapter, AMP::Mesh::Vertex, "Residual" );
 
     siloWriter->writeFile( input_file , 0 );
#endif

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



