
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
#include "vectors/Variable.h"

#include "utils/Writer.h"
#include "vectors/Vector.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/VolumeIntegralOperator.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/ElementOperationFactory.h"

#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscKrylovSolver.h"

#include "operators/NeutronicsRhs.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "operators/boundary/DirichletMatrixCorrection.h"
#include "ampmesh/Mesh.h"
#include "vectors/VectorBuilder.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "solvers/trilinos/TrilinosMLSolver.h"


void linearThermalTest(AMP::UnitTest *ut )
{
  // Input and output file names
  //  #include <string>
  std::string exeName("testTrilinosMLSolver-LinearThermalOperator-cylinder");
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
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

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
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> unusedModel;
  boost::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator 
    = boost::dynamic_pointer_cast<AMP::Operator::NeutronicsRhs>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
													       "NeutronicsOperator",
													       input_db,
													       unusedModel));

  neutronicsOperator->setTimeStep(0);

  AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar = neutronicsOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   SpecificPowerVec = AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );

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
  AMP::LinearAlgebra::Vector::shared_ptr   PowerInWattsVec = AMP::LinearAlgebra::createVector( nodalDofMap, PowerInWattsVar );

  // convert the vector of specific power to power for a given basis.
  sourceOperator->apply(nullVec, SpecificPowerVec, PowerInWattsVec, 1., 0.);

  ////////////////////////////////////
  //   CREATE THE THERMAL OPERATOR  //
  ////////////////////////////////////
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> diffusionOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																						       "DiffusionBVPOperator",
																						       input_db,
																						       transportModel));


  AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec = AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getInputVariable()  );
  AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec       = AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable() );
  AMP::LinearAlgebra::Vector::shared_ptr ResidualVec            = AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable() );

  RightHandSideVec->copyVector(PowerInWattsVec);

  boost::shared_ptr<AMP::Operator::BoundaryOperator> boundaryOp;
  boundaryOp = diffusionOperator->getBoundaryOperator();

  boundaryOp->addRHScorrection(RightHandSideVec);
  boundaryOp->setRHScorrection(RightHandSideVec);
/*
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

  double rhsNorm = RightHandSideVec->L2Norm();
  std::cout<<"RHS Norm: "<<rhsNorm<<std::endl;

  // Create the ML Solver
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver>         mlSolver(new AMP::Solver::TrilinosMLSolver(mlSolverParams));

  // Use a random initial guess?
  mlSolver->setZeroInitialGuess(false);

  // Solve the prblem.
  mlSolver->solve(RightHandSideVec, TemperatureInKelvinVec);
*/
//----------------------------------------------------------------------------------------------------
  boost::shared_ptr<AMP::Database>                 linearSolver_db = input_db->getDatabase("LinearSolver"); 

  // ---- first initialize the preconditioner
  boost::shared_ptr<AMP::Database> pcSolver_db = linearSolver_db->getDatabase("Preconditioner"); 
  boost::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> pcSolverParams(new AMP::Solver::TrilinosMLSolverParameters(pcSolver_db));
  pcSolverParams->d_pOperator = diffusionOperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(new AMP::Solver::TrilinosMLSolver(pcSolverParams));

  // initialize the linear solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));
  linearSolverParams->d_pOperator = diffusionOperator;
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pPreconditioner = pcSolver;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));
//----------------------------------------------------------------------------------------------------
  linearSolver->setZeroInitialGuess(false);
  linearSolver->solve(RightHandSideVec, TemperatureInKelvinVec);
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

/*
  // check the solution
  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator iterator = meshAdapter->beginOwnedNode(); 
  AMP::Mesh::DOFMap::shared_ptr dofmap = meshAdapter->getDOFMap( diffusionOperator->getInputVariable() );

  // The analytical solution is:  T = a + b*z + c*z*z
  //   c = -power/2
  //   b = -10*power
  //   a = 300 + 150*power


  double power = 1.;
  double radiusInMeters = 0.00534;
  double a = ;
  double b = ;
  double cal, r, sol;
  bool passes = 1;

  for( ; iterator != meshAdapter->endOwnedNode(); iterator++ ) {
    cal = TemperatureInKelvinVec->getValueByGlobalID( dofmap->getGlobalID( iterator->globalID(), 0 ) );
    r = sqrt( iterator->x()*iterator->x() + iterator->y()*iterator->y() );
    sol = function of r;
    if( fabs(cal - sol) > cal*1e-3 ) {
      passes = 0;
      ITFAILS;
    }
  }
  if( passes ) ut.passes("The linear thermal solve is verified.");
*/ 

 // Plot the results
#ifdef USE_EXT_SILO
     AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
     siloWriter->registerMesh( meshAdapter );

     siloWriter->registerVector( SpecificPowerVec,       meshAdapter, AMP::Mesh::Volume, "SpecificPower" );
     siloWriter->registerVector( PowerInWattsVec,        meshAdapter, AMP::Mesh::Vertex, "PowerInWatts" );
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


