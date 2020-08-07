
#include "AMP/materials/Material.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include <memory>
#include <string>

#include "AMP/operators/ElementOperationFactory.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/Vector.h"

#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"

#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/OperatorBuilder.h"

#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/vectors/VectorBuilder.h"

#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"


static void linearThermalTest( AMP::UnitTest *ut )
{
    // Input and output file names
    //  #include <string>
    std::string exeName( "testTrilinosMLSolver-LinearThermalOperator-cylinder" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    ////////////////////////////////////
    //    INITIALIZE THE PROBLEM      //
    ////////////////////////////////////

    // Construct a smart pointer to a new database.
    //  #include <memory>
    //  #include "AMP/utils/Database.h"


    // Fill the database from the input file.
    //  #include "AMP/utils/Database.h"
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );


    // Print from all cores into the output files
    //   #include "AMP/utils/PIO.h"
    AMP::PIO::logAllNodes( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    //--------------------------------------------------
    //   Create the Mesh.
    //--------------------------------------------------
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    std::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );
    //--------------------------------------------------

    //--------------------------------------------------
    // Create a DOF manager for a nodal vector
    //--------------------------------------------------
    int DOFsPerNode          = 1;
    int DOFsPerElement       = 8;
    int nodalGhostWidth      = 1;
    int gaussPointGhostWidth = 1;
    bool split               = true;
    AMP::Discretization::DOFManager::shared_ptr nodalDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
    AMP::Discretization::DOFManager::shared_ptr gaussPointDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Volume, gaussPointGhostWidth, DOFsPerElement, split );
    //--------------------------------------------------

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    ////////////////////////////////////
    //  CREATE THE NEUTRONICS SOURCE  //
    ////////////////////////////////////
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> unusedModel;
    std::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator =
        std::dynamic_pointer_cast<AMP::Operator::NeutronicsRhs>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "NeutronicsOperator", input_db, unusedModel ) );

    neutronicsOperator->setTimeStep( 0 );

    AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar =
        neutronicsOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr SpecificPowerVec =
        AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );

    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    /////////////////////////////////////////////////////
    //  Integrate Nuclear Source over Desnity * GeomType::Volume //
    /////////////////////////////////////////////////////

    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
    std::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator =
        std::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "VolumeIntegralOperator", input_db, stransportModel ) );

    // Create the power (heat source) vector.
    AMP::LinearAlgebra::Variable::shared_ptr PowerInWattsVar = sourceOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr PowerInWattsVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, PowerInWattsVar );

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, PowerInWattsVec );

    ////////////////////////////////////
    //   CREATE THE THERMAL OPERATOR  //
    ////////////////////////////////////
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> diffusionOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "DiffusionBVPOperator", input_db, transportModel ) );


    AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getInputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable() );

    RightHandSideVec->copyVector( PowerInWattsVec );

    std::shared_ptr<AMP::Operator::BoundaryOperator> boundaryOp;
    boundaryOp = diffusionOperator->getBoundaryOperator();

    boundaryOp->addRHScorrection( RightHandSideVec );
    boundaryOp->setRHScorrection( RightHandSideVec );

    /*
      // make sure the database on theinput file exists for the linear solver
      AMP_INSIST(input_db->keyExists("LinearSolver"),   "Key ''LinearSolver'' is missing!");

      // Read the input file onto a database.
      std::shared_ptr<AMP::Database>                 mlSolver_db   =
      input_db->getDatabase("LinearSolver");

      // Fill in the parameters fo the class with the info on the database.
      std::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams (new
      AMP::Solver::SolverStrategyParameters(mlSolver_db));

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
      std::shared_ptr<AMP::Solver::TrilinosMLSolver>         mlSolver(new
      AMP::Solver::TrilinosMLSolver(mlSolverParams));

      // Use a random initial guess?
      mlSolver->setZeroInitialGuess(false);

      // Solve the prblem.
      mlSolver->solve(RightHandSideVec, TemperatureInKelvinVec);
    */
    //----------------------------------------------------------------------------------------------------
    std::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase( "LinearSolver" );

    // ---- first initialize the preconditioner
    std::shared_ptr<AMP::Database> pcSolver_db = linearSolver_db->getDatabase( "Preconditioner" );
    std::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> pcSolverParams(
        new AMP::Solver::TrilinosMLSolverParameters( pcSolver_db ) );
    pcSolverParams->d_pOperator = diffusionOperator;
    std::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(
        new AMP::Solver::TrilinosMLSolver( pcSolverParams ) );

    // initialize the linear solver
    std::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(
        new AMP::Solver::PetscKrylovSolverParameters( linearSolver_db ) );
    linearSolverParams->d_pOperator       = diffusionOperator;
    linearSolverParams->d_comm            = globalComm;
    linearSolverParams->d_pPreconditioner = pcSolver;
    std::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(
        new AMP::Solver::PetscKrylovSolver( linearSolverParams ) );
    //----------------------------------------------------------------------------------------------------
    linearSolver->setZeroInitialGuess( false );
    linearSolver->solve( RightHandSideVec, TemperatureInKelvinVec );
    // Compute the residual
    diffusionOperator->residual( RightHandSideVec, TemperatureInKelvinVec, ResidualVec );

    // Check the L2 norm of the final residual.
    double finalResidualNorm = ResidualVec->L2Norm(ResidualVec);
    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    if ( finalResidualNorm > 10.0 ) {
        ut->failure( "TrilinosMLSolver successfully solves a linear thermal problem with a nuclear "
                     "source term." );
    } else {
        ut->passes( "TrilinosMLSolver successfully solves a linear thermal problem with a nuclear "
                    "source term." );
    }

/*
  // check the solution
  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator iterator = meshAdapter->beginOwnedNode();
  AMP::Mesh::DOFMap::shared_ptr dofmap = meshAdapter->getDOFMap(
  diffusionOperator->getInputVariable() );

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
    cal = TemperatureInKelvinVec->getValueByGlobalID( dofmap->getGlobalID( iterator->globalID(), 0 )
  );
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
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( meshAdapter );

    siloWriter->registerVector(
        SpecificPowerVec, meshAdapter, AMP::Mesh::GeomType::Volume, "SpecificPower" );
    siloWriter->registerVector(
        PowerInWattsVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "PowerInWatts" );
    siloWriter->registerVector(
        TemperatureInKelvinVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "TemperatureInKelvin" );
    siloWriter->registerVector( ResidualVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Residual" );

    siloWriter->writeFile( input_file, 0 );
#endif

    input_db.reset();

    ut->passes( exeName );
}


int testTrilinosMLSolver_LinearThermalOperator_cylinder( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    linearThermalTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
