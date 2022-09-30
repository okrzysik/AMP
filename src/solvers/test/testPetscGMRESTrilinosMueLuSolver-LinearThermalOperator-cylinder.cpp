#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"

#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
#include "AMP/operators/ElementOperationFactory.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/solvers/trilinos/muelu/TrilinosMueLuSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <memory>


void linearThermalTest( AMP::UnitTest *ut )
{
    // Input and output file names
    std::string exeName( "testPetscGMRESTrilinosMueLuSolver-LinearThermalOperator-cylinder" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    // Fill the database from the input file.
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );


    // Print from all cores into the output files
    AMP::logAllNodes( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Create the Mesh
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
    //  CREATE THE NEUTRONICS SOURCE  //
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> unusedModel;
    auto neutronicsOperator = std::dynamic_pointer_cast<AMP::Operator::NeutronicsRhs>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "NeutronicsOperator", input_db, unusedModel ) );

    neutronicsOperator->setTimeStep( 0 );

    auto SpecificPowerVar = neutronicsOperator->getOutputVariable();
    auto SpecificPowerVec = AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );

    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    //  Integrate Nuclear Source over Desnity * Volume //
    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
    auto sourceOperator = std::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "VolumeIntegralOperator", input_db, stransportModel ) );

    // Create the power (heat source) vector.
    auto PowerInWattsVar = sourceOperator->getOutputVariable();
    auto PowerInWattsVec = AMP::LinearAlgebra::createVector( nodalDofMap, PowerInWattsVar );

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, PowerInWattsVec );

    //   CREATE THE THERMAL OPERATOR  //
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

    RightHandSideVec->copyVector( PowerInWattsVec );

    auto boundaryOp = diffusionOperator->getBoundaryOperator();

    boundaryOp->addRHScorrection( RightHandSideVec );
    boundaryOp->setRHScorrection( RightHandSideVec );

    auto linearSolver_db = input_db->getDatabase( "LinearSolver" );

    // ---- first initialize the preconditioner
    auto pcSolver_db = linearSolver_db->getDatabase( "Preconditioner" );
    auto pcSolverParams =
        std::make_shared<AMP::Solver::TrilinosMueLuSolverParameters>( pcSolver_db );
    pcSolverParams->d_pOperator = diffusionOperator;
    auto pcSolver = std::make_shared<AMP::Solver::TrilinosMueLuSolver>( pcSolverParams );

    // initialize the linear solver
    auto linearSolverParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( linearSolver_db );
    linearSolverParams->d_pOperator       = diffusionOperator;
    linearSolverParams->d_comm            = globalComm;
    linearSolverParams->d_pPreconditioner = pcSolver;
    auto linearSolver = std::make_shared<AMP::Solver::PetscKrylovSolver>( linearSolverParams );
    //----------------------------------------------------------------------------------------------------
    linearSolver->setZeroInitialGuess( false );
    linearSolver->apply( RightHandSideVec, TemperatureInKelvinVec );
    // Compute the residual
    diffusionOperator->residual( RightHandSideVec, TemperatureInKelvinVec, ResidualVec );

    // Check the L2 norm of the final residual.
    double finalResidualNorm = static_cast<double>( ResidualVec->L2Norm() );
    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    if ( finalResidualNorm > 10.0 ) {
        ut->failure( "TrilinosMueLuSolver could NOT successfully solve a linear thermal problem "
                     "with a nuclear "
                     "source term." );
    } else {
        ut->passes(
            "TrilinosMueLuSolver successfully solves a linear thermal problem with a nuclear "
            "source term." );
    }

    // Plot the results
    auto siloWriter = AMP::IO::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( meshAdapter );
    siloWriter->registerVector(
        SpecificPowerVec, meshAdapter, AMP::Mesh::GeomType::Cell, "SpecificPower" );
    siloWriter->registerVector(
        PowerInWattsVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "PowerInWatts" );
    siloWriter->registerVector(
        TemperatureInKelvinVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "TemperatureInKelvin" );
    siloWriter->registerVector( ResidualVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Residual" );
    siloWriter->writeFile( input_file, 0 );

    input_db.reset();

    ut->passes( exeName );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    linearThermalTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
