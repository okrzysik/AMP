#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshParameters.h"

#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
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
#include "AMP/solvers/trilinos/muelu/TrilinosMueLuSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <memory>
#include <string>


void linearThermalTest( AMP::UnitTest *ut, std::string exeName )
{
    double t1;
    // Input and output file names
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    size_t N_error0        = ut->NumFailLocal();

    // Fill the database from the input file.
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Print from all cores into the output files
    AMP::logAllNodes( log_file );

    // Create the Mesh
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db   = input_db->getDatabase( "Mesh" );
    auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );

    // Create a DOF manager for a nodal vector
    int DOFsPerNode          = 1;
    int DOFsPerElement       = 8;
    int nodalGhostWidth      = 1;
    int gaussPointGhostWidth = 1;
    bool split               = true;
    auto nodalDofMap         = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
    auto gaussPointDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Volume, gaussPointGhostWidth, DOFsPerElement, split );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    //  CREATE THE NEUTRONICS SOURCE
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    auto neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    auto neutronicsParams =
        std::make_shared<AMP::Operator::NeutronicsRhsParameters>( neutronicsOp_db );
    auto neutronicsOperator = std::make_shared<AMP::Operator::NeutronicsRhs>( neutronicsParams );

    auto SpecificPowerVar = neutronicsOperator->getOutputVariable();
    auto SpecificPowerVec = AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );
    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    //  Integrate Nuclear Rhs over Desnity * Volume

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
    RightHandSideVec->setToScalar( 0.0 );

    //   Add the boundary conditions corrections
    RightHandSideVec->copyVector( PowerInWattsVec );
    diffusionOperator->modifyRHSvector( RightHandSideVec );
    std::cout << "RHS Norm 1: " << RightHandSideVec->L2Norm() << std::endl;
    std::cout << "RHS Norm 2: " << PowerInWattsVec->L2Norm() << std::endl;


    /////////////////////////////////////////////
    // make sure the database on theinput file exists for the linear solver
    AMP_INSIST( input_db->keyExists( "LinearSolver" ), "Key ''LinearSolver'' is missing!" );

    // Read the input file onto a database.
    std::shared_ptr<AMP::Database> mlSolver_db = input_db->getDatabase( "LinearSolver" );

    // Fill in the parameters fo the class with the info on the database.
    auto mlSolverParams = std::make_shared<AMP::Solver::SolverStrategyParameters>( mlSolver_db );

    // Define the operature to be used by the Solver.
    mlSolverParams->d_pOperator = diffusionOperator;

    // Set initial guess
    TemperatureInKelvinVec->setToScalar( 1.0 );

    // Check the initial L2 norm of the solution
    std::cout << "Initial Solution Norm: " << TemperatureInKelvinVec->L2Norm() << std::endl;
    std::cout << "RHS Norm: " << RightHandSideVec->L2Norm() << std::endl;

    // Create the ML Solver
    auto mlSolver = std::make_shared<AMP::Solver::TrilinosMueLuSolver>( mlSolverParams );

    // Use a random initial guess?
    mlSolver->setZeroInitialGuess( false );

    AMP::pout << "RHS Max: " << RightHandSideVec->max() << std::endl;
    AMP::pout << "RHS Min: " << RightHandSideVec->min() << std::endl;
    AMP::pout << "RHS L1-norm: " << RightHandSideVec->L1Norm() << std::endl;
    AMP::pout << "RHS L2-norm: " << RightHandSideVec->L2Norm() << std::endl;

    // Solve the prblem.
    mlSolver->apply( RightHandSideVec, TemperatureInKelvinVec );

    AMP::pout << "Solution Max: " << TemperatureInKelvinVec->max() << std::endl;
    AMP::pout << "Solution Min: " << TemperatureInKelvinVec->min() << std::endl;
    AMP::pout << "Solution L1-norm: " << TemperatureInKelvinVec->L1Norm() << std::endl;
    AMP::pout << "Solution L2-norm: " << TemperatureInKelvinVec->L2Norm() << std::endl;

    // Compute the residual
    diffusionOperator->residual( RightHandSideVec, TemperatureInKelvinVec, ResidualVec );

    // Check the L2 norm of the final residual.
    double finalResidualNorm = static_cast<double>( ResidualVec->L2Norm() );
    AMP::pout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    if ( finalResidualNorm > 10.0 ) {
        ut->failure( "TrilinosMueLuSolver could not solve a linear thermal problem succesfully" );
    }

    // Plot the results
    auto globalComm = AMP::AMP_MPI( AMP_COMM_WORLD );
#ifdef USE_EXT_SILO
    auto siloWriter = AMP::IO::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( meshAdapter );
    siloWriter->registerVector(
        TemperatureInKelvinVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "TemperatureInKelvin" );
    siloWriter->registerVector( ResidualVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Residual" );

    siloWriter->writeFile( input_file, 0 );
#endif

    input_db.reset();

    if ( N_error0 == ut->NumFailLocal() )
        ut->passes( exeName );
    else
        ut->failure( exeName );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;

    exeNames.push_back( "testTrilinosMueLuSolver-LinearThermalOperator-2_HALDEN" );
    exeNames.push_back( "testTrilinosMueLuSolver-LinearThermalOperator-cylinder" );
    exeNames.push_back( "testTrilinosMueLuSolver-LinearThermalOperator-shell" );
    //    exeNames.push_back( "testTrilinosMueLuSolver-LinearThermalOperator-2_HALDEN_clad" );

    for ( auto &exeName : exeNames )
        linearThermalTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
