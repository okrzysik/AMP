#include "materials/Material.h"
#include "operators/NeutronicsRhs.h"
#include "utils/AMPManager.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/Utilities.h"
#include "utils/shared_ptr.h"
#include "vectors/Variable.h"
#include <string>

#include "operators/ElementOperationFactory.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/libmesh/VolumeIntegralOperator.h"
#include "utils/Writer.h"
#include "vectors/Vector.h"

#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/boundary/libmesh/NeumannVectorCorrection.h"
#include "operators/boundary/libmesh/RobinMatrixCorrection.h"

#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "solvers/SolverFactory.h"
#include "solvers/KrylovSolverParameters.h"

AMP::shared_ptr<AMP::Solver::SolverStrategy>
buildSolver ( const AMP::shared_ptr<AMP::InputDatabase> &input_db,
              const std::string &solver_name,                            
              const AMP::AMP_MPI &comm,
              AMP::shared_ptr<AMP::Operator::Operator> &op )
{

    AMP::shared_ptr<AMP::Solver::SolverStrategy> solver;
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> parameters;
    
    AMP_INSIST( input_db->keyExists( solver_name ), "Key " + solver_name + " is missing!" );
    
    const auto &db = input_db->getDatabase( solver_name );
    
    if ( db->keyExists( "name" ) ) {
        
        auto name = db->getString( "name" );

        if ( ( name == "GMRESSolver" )
          || ( name == "CGSolver" )
          || ( name == "BiCGSTABSolver" )
          || ( name == "TFQMRSolver" ) ) {

            // check if we need to construct a preconditioner
            auto use_preconditioner = db->getBoolWithDefault( "use_preconditioner", false );
            AMP::shared_ptr<AMP::Solver::SolverStrategy> pcSolver;

            if ( use_preconditioner ) {

                auto pc_name = db->getStringWithDefault( "pc_name", "Preconditioner" );

                pcSolver = buildSolver( input_db,
                                        pc_name,
                                        comm,
                                        op );
                
                AMP_INSIST( pcSolver.get() != nullptr, "null preconditioner"); 
            }
            
            auto params = AMP::make_shared<AMP::Solver::KrylovSolverParameters> ( db );
            params->d_comm = comm;
            params->d_pPreconditioner = pcSolver;
            parameters = params;
            
        } else {
            parameters = AMP::make_shared<AMP::Solver::SolverStrategyParameters> ( db );
        }
        
        AMP_INSIST( parameters != nullptr, "null parameter object" );
        parameters->d_pOperator = op;

    } else {
        AMP_ERROR("Key name does not exist in solver database");
    }

    solver = AMP::Solver::SolverFactory::create( parameters );

    return solver;

}

void linearThermalTest( AMP::UnitTest *ut, std::string inputFileName )
{
    // double t1;
    // Input and output file names
    //  #include <string>
    std::string input_file = inputFileName;
    std::string log_file   = "output_" + inputFileName;
    ////////////////////////////////////
    //    INITIALIZE THE PROBLEM      //
    ////////////////////////////////////

    size_t N_error0 = ut->NumFailLocal();

    // Construct a smart pointer to a new database.
    //  #include "utils/shared_ptr.h"
    //  #include "utils/InputDatabase.h"
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );

    // Fill the database from the input file.
    //  #include "utils/InputManager.h"
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );


    // Print from all cores into the output files
    //   #include "utils/PIO.h"
    AMP::PIO::logAllNodes( log_file );

    //--------------------------------------------------
    //   Create the Mesh.
    //--------------------------------------------------
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );
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
    AMP::shared_ptr<AMP::Database> neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    AMP::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(
        new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ) );
    AMP::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(
        new AMP::Operator::NeutronicsRhs( neutronicsParams ) );

    AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar =
        neutronicsOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr SpecificPowerVec =
        AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );

    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    /////////////////////////////////////////////////////
    //  Integrate Nuclear Source over Desnity * Volume //
    /////////////////////////////////////////////////////

    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "VolumeIntegralOperator", input_db, stransportModel ) );

    // Create the power (heat source) vector.
    AMP::LinearAlgebra::Variable::shared_ptr PowerInWattsVar = sourceOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr PowerInWattsVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, PowerInWattsVar );
    PowerInWattsVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, PowerInWattsVec );

    ////////////////////////////////////////
    //   CREATE THE THERMAL BVP OPERATOR  //
    ////////////////////////////////////////
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
    auto linearOperator = AMP::Operator::OperatorBuilder::createOperator( meshAdapter,
                                                                             "DiffusionBVPOperator",
                                                                             input_db,
                                                                             transportModel );
    
    auto diffusionOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>( linearOperator );

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getInputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable() );

    RightHandSideVec->setToScalar( 0.0 );
    double rhsNorm = RightHandSideVec->L2Norm();

    ///////////////////////////////////////////////
    //   Add the boundary conditions corrections //
    ///////////////////////////////////////////////

    AMP::LinearAlgebra::Vector::shared_ptr boundaryOpCorrectionVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable() );

    AMP::Operator::Operator::shared_ptr boundaryOp;
    boundaryOp = diffusionOperator->getBoundaryOperator();

    ( AMP::dynamic_pointer_cast<AMP::Operator::BoundaryOperator>( boundaryOp ) )
        ->addRHScorrection( boundaryOpCorrectionVec );

    RightHandSideVec->subtract( PowerInWattsVec, boundaryOpCorrectionVec );

    rhsNorm = RightHandSideVec->L2Norm();
    std::cout << "RHS Norm after BC Correction " << rhsNorm << std::endl;

    rhsNorm = RightHandSideVec->L2Norm();
    std::cout << "RHS Norm 1: " << rhsNorm << std::endl;
    rhsNorm = PowerInWattsVec->L2Norm();
    std::cout << "RHS Norm 2: " << rhsNorm << std::endl;
    rhsNorm = boundaryOpCorrectionVec->L2Norm();
    std::cout << "RHS Norm 3: " << rhsNorm << std::endl;

    /////////////////////////////////////////////
    // make sure the database on theinput file exists for the linear solver
    AMP_INSIST( input_db->keyExists( "LinearSolver" ), "Key ''LinearSolver'' is missing!" );

    auto comm = AMP::AMP_MPI( AMP_COMM_WORLD );
    auto linearSolver = buildSolver(  input_db, "LinearSolver", comm, linearOperator );

    // Set initial guess
    TemperatureInKelvinVec->setToScalar( 1.0 );

    // Check the initial L2 norm of the solution
    double initSolNorm = TemperatureInKelvinVec->L2Norm();
    std::cout << "Initial Solution Norm: " << initSolNorm << std::endl;

    rhsNorm = RightHandSideVec->L2Norm();
    std::cout << "RHS Norm: " << rhsNorm << std::endl;

    // Use a random initial guess?
    linearSolver->setZeroInitialGuess( false );

    // Solve the prblem.
    linearSolver->solve( RightHandSideVec, TemperatureInKelvinVec );

    // Compute the residual
    diffusionOperator->residual( RightHandSideVec, TemperatureInKelvinVec, ResidualVec );

    // Check the L2 norm of the final residual.
    double finalResidualNorm = ResidualVec->L2Norm();
    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    if ( finalResidualNorm > 10.0 ) {
        ut->failure( "TrilinosMueLuSolver does not solve a linear thermal problem with a nuclear "
                     "source term." );
    } 

    // Plot the results
    AMP::AMP_MPI globalComm = AMP::AMP_MPI( AMP_COMM_WORLD );
#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( meshAdapter );

    siloWriter->registerVector(
        TemperatureInKelvinVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "TemperatureInKelvin" );
    siloWriter->registerVector( ResidualVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Residual" );

    siloWriter->writeFile( input_file, 0 );
#endif

    input_db.reset();

    if ( N_error0 == ut->NumFailLocal() )
        ut->passes( inputFileName );
    else
        ut->failure( inputFileName );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    AMP::Solver::registerSolverFactories();

    std::vector<std::string> files;
    
    if ( argc > 1 ) {

        files.push_back( argv[1] );

    } else {
    
        files.push_back( "input_testLinearSolvers-LinearThermalRobin-GMRES" );
        files.push_back( "input_testLinearSolvers-LinearThermalRobin-BiCGSTAB" );
        
#ifdef USE_EXT_HYPRE
        files.push_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG" );
        files.push_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRES" );
        files.push_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-BiCGSTAB" );
#endif

#ifdef USE_TRILINOS_ML
        files.push_back( "input_testLinearSolvers-LinearThermalRobin-ML" );
        files.push_back( "input_testLinearSolvers-LinearThermalRobin-ML-GMRES" );
        files.push_back( "input_testLinearSolvers-LinearThermalRobin-ML-BiCGSTAB" );
#endif

#ifdef USE_TRILINOS_MUELU
        files.push_back( "input_testLinearSolvers-LinearThermalRobin-MueLu" );
        files.push_back( "input_testLinearSolvers-LinearThermalRobin-MueLu-GMRES" );
        files.push_back( "input_testLinearSolvers-LinearThermalRobin-MueLu-BiCGSTAB" );
#endif

    }

    for ( auto &file : files )
        linearThermalTest( &ut, file );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
