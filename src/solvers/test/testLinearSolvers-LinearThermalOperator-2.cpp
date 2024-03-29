#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/operators/ElementOperationFactory.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/libmesh/NeumannVectorCorrection.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/solvers/SolverFactory.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/solvers/testHelpers/SolverTestParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <memory>

std::shared_ptr<AMP::Solver::SolverStrategy>
buildSolver( std::shared_ptr<AMP::Database> input_db,
             const std::string &solver_name,
             const AMP::AMP_MPI &comm,
             std::shared_ptr<AMP::Operator::Operator> op )
{

    AMP_INSIST( input_db->keyExists( solver_name ), "Key " + solver_name + " is missing!" );

    auto db = input_db->getDatabase( solver_name );
    AMP_INSIST( db->keyExists( "name" ), "Key name does not exist in solver database" );

    auto parameters         = std::make_shared<AMP::Solver::SolverStrategyParameters>( db );
    parameters->d_pOperator = op;
    parameters->d_comm      = comm;
    parameters->d_global_db = input_db;

    return AMP::Solver::SolverFactory::create( parameters );
}

void linearThermalTest( AMP::UnitTest *ut,
                        const std::string &input_file,
                        std::shared_ptr<AMP::Database> input_db,
                        std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
                        std::shared_ptr<AMP::LinearAlgebra::Vector> powerVec )
{

    AMP_ASSERT( input_db && meshAdapter && powerVec );

    auto linearSolverName = input_db->getDatabase( "LinearSolver" )->getString( "name" );

    std::string pcName;
    if ( input_db->keyExists( "Preconditioner" ) )
        pcName = input_db->getDatabase( "Preconditioner" )->getString( "name" );
    else
        pcName = "No PC";

    // Create a DOF manager for a nodal vector
    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    auto nodalDofMap    = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );

    ////////////////////////////////////
    //   CREATE THE THERMAL OPERATOR  //
    ////////////////////////////////////
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
    auto diffusionOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "DiffusionBVPOperator", input_db, transportModel ) );
    auto TemperatureInKelvinVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getInputVariable() );
    auto RightHandSideVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable() );
    auto ResidualVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable() );
    RightHandSideVec->zero();

    //   Add the boundary conditions corrections //
    RightHandSideVec->copyVector( powerVec );
    diffusionOperator->modifyRHSvector( RightHandSideVec );
    std::cout << "RHS Norm 1: " << RightHandSideVec->L2Norm() << std::endl;
    std::cout << "RHS Norm 2: " << powerVec->L2Norm() << std::endl;

    // Set initial guess
    TemperatureInKelvinVec->setToScalar( 1.0 );

    // Check the initial L2 norm of the solution
    std::cout << "Initial Solution Norm: " << TemperatureInKelvinVec->L2Norm() << std::endl;
    std::cout << "RHS Norm: " << RightHandSideVec->L2Norm() << std::endl;

    auto comm         = AMP::AMP_MPI( AMP_COMM_WORLD );
    auto linearSolver = buildSolver( input_db, "LinearSolver", comm, diffusionOperator );

    AMP::pout << "RHS Max: " << RightHandSideVec->max() << std::endl;
    AMP::pout << "RHS Min: " << RightHandSideVec->min() << std::endl;
    AMP::pout << "RHS L1-norm: " << RightHandSideVec->L1Norm() << std::endl;
    AMP::pout << "RHS L2-norm: " << RightHandSideVec->L2Norm() << std::endl;

    // Solve the problem.
    linearSolver->setZeroInitialGuess( false );
    linearSolver->apply( RightHandSideVec, TemperatureInKelvinVec );

    AMP::pout << "Solution Max: " << TemperatureInKelvinVec->max() << std::endl;
    AMP::pout << "Solution Min: " << TemperatureInKelvinVec->min() << std::endl;
    AMP::pout << "Solution L1-norm: " << TemperatureInKelvinVec->L1Norm() << std::endl;
    AMP::pout << "Solution L2-norm: " << TemperatureInKelvinVec->L2Norm() << std::endl;
    AMP::pout << "======================================================================="
              << std::endl;

    // Compute the residual
    diffusionOperator->residual( RightHandSideVec, TemperatureInKelvinVec, ResidualVec );

    // Check the L2 norm of the final residual.
    double finalResidualNorm = static_cast<double>( ResidualVec->L2Norm() );
    AMP::pout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    auto combo = linearSolverName + " with PC " + pcName;
    if ( finalResidualNorm < 10 ) {
        ut->passes( combo + " passes linear thermal problem with input " + input_file );
    } else {
        ut->passes( combo + " fails linear thermal problem with input " + input_file );
    }

#ifdef AMP_USE_SILO
    // Plot the results
    auto siloWriter = AMP::IO::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( meshAdapter );
    siloWriter->registerVector(
        TemperatureInKelvinVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "TemperatureInKelvin" );
    siloWriter->registerVector( ResidualVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Residual" );
    siloWriter->writeFile( input_file, 0 );
#endif
}

void linearThermalTest( AMP::UnitTest *ut, const std::string &inputFile )
{
    double t1;
    // Input and output file names
    std::string input_file = inputFile;
    std::string log_file   = "output_" + inputFile;

    // Fill the database from the input file
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Print from all cores into the output files
    AMP::logAllNodes( log_file );

    //   Create the Mesh.
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db   = input_db->getDatabase( "Mesh" );
    auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::MeshFactory::create( mgrParams );

    // Create a DOF manager for a nodal vector
    int DOFsPerNode          = 1;
    int DOFsPerElement       = 8;
    int nodalGhostWidth      = 1;
    int gaussPointGhostWidth = 1;
    bool split               = true;
    auto nodalDofMap         = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
    auto gaussPointDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Cell, gaussPointGhostWidth, DOFsPerElement, split );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    // CREATE THE NEUTRONICS SOURCE
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    auto neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    auto neutronicsParams =
        std::make_shared<AMP::Operator::NeutronicsRhsParameters>( neutronicsOp_db );
    auto neutronicsOperator = std::make_shared<AMP::Operator::NeutronicsRhs>( neutronicsParams );
    auto SpecificPowerVar   = neutronicsOperator->getOutputVariable();
    auto SpecificPowerVec = AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );
    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    // Integrate Nuclear Rhs over Density * Volume //
    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
    auto sourceOperator = std::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "VolumeIntegralOperator", input_db, stransportModel ) );

    // Create the power (heat source) vector.
    auto PowerInWattsVar = sourceOperator->getOutputVariable();
    auto PowerInWattsVec = AMP::LinearAlgebra::createVector( nodalDofMap, PowerInWattsVar );
    PowerInWattsVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, PowerInWattsVec );

    t1 = static_cast<double>( SpecificPowerVec->L2Norm() );
    std::cout << "n1 = " << t1 << std::endl;
    t1 = static_cast<double>( PowerInWattsVec->L2Norm() );
    std::cout << "n1 = " << t1 << std::endl;

    std::vector<std::pair<std::string, std::string>> solvers{
#ifdef AMP_USE_HYPRE
        { "CG", "BoomerAMG" },
        { "GMRES", "BoomerAMG" },
        { "FGMRES", "BoomerAMG" },
        { "BiCGSTAB", "BoomerAMG" },
        { "TFQMR", "BoomerAMG" },
        { "BoomerAMG", "NoPC" },
        { "HyprePCG", "NoPC" },
        { "HyprePCG", "BoomerAMG" },
    #ifdef AMP_USE_PETSC
        { "PetscFGMRES", "BoomerAMG" },
    #endif
#endif
#ifdef AMP_USE_TRILINOS_ML
        { "CG", "ML" },
        { "GMRES", "ML" },
        { "FGMRES", "ML" },
        { "BiCGSTAB", "ML" },
        { "TFQMR", "ML" },
    #ifdef AMP_USE_PETSC
        { "PetscFGMRES", "ML" },
    #endif
        { "ML", "NoPC" },
#endif
#ifdef AMP_USE_PETSC
        { "PetscFGMRES", "NoPC" },
#endif
        { "CG", "NoPC" },
        { "GMRES", "NoPC" },
        { "FGMRES", "NoPC" },
        { "BiCGSTAB", "NoPC" },
        { "TFQMR", "NoPC" }
    };

    for ( auto &[primary, nested] : solvers ) {
        std::shared_ptr<AMP::Database> db = input_db->cloneDatabase();
        auto use_nested                   = ( nested == "NoPC" ) ? false : true;
        db->putDatabase(
            "LinearSolver",
            AMP::Solver::Test::SolverParameters::getParameters( primary, use_nested ) );
        if ( use_nested ) {
            db->putDatabase(
                "Preconditioner",
                AMP::Solver::Test::SolverParameters::getParameters( nested, use_nested ) );
        }

        std::string banner;
        if ( use_nested )
            banner = "Running " + primary + " with PC " + nested + " on " + input_file;
        else
            banner = "Running " + primary + " on " + input_file;
        AMP::pout << banner << std::endl;
        linearThermalTest( ut, input_file + primary + nested, db, meshAdapter, PowerInWattsVec );
    }
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    AMP::Solver::registerSolverFactories();

    std::vector<std::string> inputFiles;

    if ( argc > 1 ) {

        inputFiles.emplace_back( argv[1] );

    } else {

        inputFiles.emplace_back( "input_LinearThermalOperator-2_HALDEN" );
        inputFiles.emplace_back( "input_LinearThermalOperator-2-cylinder" );
        inputFiles.emplace_back( "input_LinearThermalOperator-2-shell" );
        //    inputFiles.push_back( "input_testBoomerAMGSolver-LinearThermalOperator-2_HALDEN_clad"
        //    );
    }

    for ( auto &inputFile : inputFiles )
        linearThermalTest( &ut, inputFile );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
