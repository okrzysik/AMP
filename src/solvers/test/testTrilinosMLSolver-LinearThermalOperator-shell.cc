#include "AMP/materials/Material.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/Variable.h"
#include <string>

#include "AMP/operators/DirichletMatrixCorrection.h"
#include "AMP/operators/ElementOperationFactory.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/Vector.h"

#include "../TrilinosMLSolver.h"


void linearThermalTest( AMP::UnitTest *ut )
{
    // Input and output file names
    //  #include <string>
    std::string exeName( "testTrilinosMLSolver-LinearThermalOperator-shell" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    ////////////////////////////////////
    //    INITIALIZE THE PROBLEM      //
    ////////////////////////////////////
    // Startup AMP by initializing MPI, IO, and exception handlers.
    //  #include "AMP/utils/AMPManager.h"
    AMP::AMPManager::startup();

    // Create the map to get an available material from a string.
    //  #include "AMP/materials/Material.h"
    AMP::Materials::initialize();

    // Construct a smart pointer to a new database.
    //  #include "AMP/utils/shared_ptr.h"
    //  #include "AMP/utils/InputDatabase.h"
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );

    // Fill the database from the input file.
    //  #include "AMP/utils/InputManager.h"
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );


    // Print from all cores into the output files
    //   #include "AMP/utils/PIO.h"
    AMP::PIO::logAllNodes( log_file );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );

    // Construct a mesh manager which reads in the fuel mesh
    AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams(
        new AMP::Mesh::MeshManagerParameters( input_db ) );
    AMP::Mesh::MeshManager::shared_ptr manager( new AMP::Mesh::MeshManager( mgrParams ) );
    AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh( "clad" );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    ////////////////////////////////////
    //  CREATE THE NEUTRONICS SOURCE  //
    ////////////////////////////////////
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    AMP::shared_ptr<AMP::Database> neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> unusedModel;
    AMP::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NeutronicsRhs>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "NeutronicsOperator", input_db, unusedModel ) );

    neutronicsOperator->setTimeStep( 0 );

    AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar =
        neutronicsOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr SpecificPowerVec =
        meshAdapter->createVector( SpecificPowerVar );

    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    /////////////////////////////////////////////////////
    //  Integrate Nuclear Source over Desnity * GeomType::Volume //
    /////////////////////////////////////////////////////

    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
    AMP::shared_ptr<AMP::Database> sourceDatabase =
        input_db->getDatabase( "VolumeIntegralOperator" );
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "VolumeIntegralOperator", input_db, stransportModel ) );

    // Create the power (heat source) vector.
    AMP::LinearAlgebra::Variable::shared_ptr PowerInWattsVar = sourceOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr PowerInWattsVec =
        meshAdapter->createVector( PowerInWattsVar );

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, PowerInWattsVec );

    ////////////////////////////////////
    //   CREATE THE THERMAL OPERATOR  //
    ////////////////////////////////////
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
    AMP::shared_ptr<AMP::InputDatabase> bvpDatabase = AMP::dynamic_pointer_cast<AMP::InputDatabase>(
        input_db->getDatabase( "DiffusionBVPOperator" ) );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> diffusionOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "DiffusionBVPOperator", input_db, transportModel ) );


    AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec =
        meshAdapter->createVector( diffusionOperator->getInputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec =
        meshAdapter->createVector( diffusionOperator->getOutputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec =
        meshAdapter->createVector( diffusionOperator->getOutputVariable() );

    RightHandSideVec->copyVector( PowerInWattsVec );

    AMP::shared_ptr<AMP::Operator::BoundaryOperator> boundaryOp;
    boundaryOp = diffusionOperator->getBoundaryOperator();

    boundaryOp->addRHScorrection( RightHandSideVec );
    boundaryOp->setRHScorrection( RightHandSideVec );

    // make sure the database on theinput file exists for the linear solver
    AMP_INSIST( input_db->keyExists( "LinearSolver" ), "Key ''LinearSolver'' is missing!" );

    // Read the input file onto a database.
    AMP::shared_ptr<AMP::Database> mlSolver_db = input_db->getDatabase( "LinearSolver" );

    // Fill in the parameters fo the class with the info on the database.
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams(
        new AMP::Solver::SolverStrategyParameters( mlSolver_db ) );

    // Define the operature to be used by the Solver.
    mlSolverParams->d_pOperator = diffusionOperator;

    // Set initial guess
    TemperatureInKelvinVec->setToScalar( 1.0 );

    // Check the initial L2 norm of the solution
    double initSolNorm = TemperatureInKelvinVec->L2Norm();
    std::cout << "Initial Solution Norm: " << initSolNorm << std::endl;

    double rhsNorm = RightHandSideVec->L2Norm();
    std::cout << "RHS Norm: " << rhsNorm << std::endl;

    // Create the ML Solver
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> mlSolver(
        new AMP::Solver::TrilinosMLSolver( mlSolverParams ) );

    // Use a random initial guess?
    mlSolver->setZeroInitialGuess( false );

    // Solve the prblem.
    mlSolver->solve( RightHandSideVec, TemperatureInKelvinVec );

    // Compute the residual
    diffusionOperator->residual( RightHandSideVec, TemperatureInKelvinVec, ResidualVec );

    // Check the L2 norm of the final residual.
    double finalResidualNorm = ResidualVec->L2Norm();
    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    if ( finalResidualNorm > 10.0 ) {
        ITFAILS;
    } else {
        ut.passes( "TrilinosMLSolver successfully solves a linear thermal problem with a nuclear "
                   "source term." );
    }

    // check the solution
    AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator iterator = meshAdapter->beginOwnedNode();
    AMP::Mesh::DOFMap::shared_ptr dofmap =
        meshAdapter->getDOFMap( diffusionOperator->getInputVariable() );

    // The analytical solution is:  T = T[r1] + (T[r2]-T[r1])*ln(r/r1)/a + b*[ (r*r-r1*r1) -
    // c/a*ln(r/r1) ]
    //   c = r2*r2 - r1*r1
    //   b = power/(4 k)
    //   a = ln(r/r1) / ln(r2/r1) ... assume r1 < r2

    double power           = 10000000.0;
    double radius1InMeters = 0.00545465;
    double radius2InMeters = 0.00639445;
    double T1InK           = 500.;
    double T2InK           = 300.;
    double a               = log( radius2InMeters / radius1InMeters );
    double b               = power / 4.;
    double c               = radius2InMeters * radius2InMeters - radius1InMeters * radius1InMeters;
    bool passes            = 1;
    double cal, r, sol;

    for ( ; iterator != meshAdapter->endOwnedNode(); iterator++ ) {
        cal = TemperatureInKelvinVec->getValueByGlobalID(
            dofmap->getGlobalID( iterator->globalID(), 0 ) );
        r   = sqrt( iterator->x() * iterator->x() + iterator->y() * iterator->y() );
        sol = T1InK + ( T2InK - T1InK ) / a * log( r / radius1InMeters ) +
              b * ( ( r * r - radius1InMeters * radius1InMeters ) -
                    c / a * log( r / radius1InMeters ) );
        if ( fabs( cal - sol ) > cal * 1e-3 ) {
            passes = 0;
            ITFAILS;
        }
    }
    if ( passes )
        ut.passes( "The linear thermal solve is verified." );

    // Plot the results
    if ( AMP::AMP_MPI::getNodes() == 1 ) {
#ifdef USE_EXT_SILO
        AMP::LinearAlgebra::Variable::shared_ptr tmpVar1 = PowerInWattsVec->getVariable();
        tmpVar1->setName( "PowerInWatts" );
        meshAdapter->registerVectorAsData( PowerInWattsVec );

        tmpVar1 = SpecificPowerVec->getVariable();
        tmpVar1->setName( "SpecificPower" );
        meshAdapter->registerVectorAsData( SpecificPowerVec );

        tmpVar1 = TemperatureInKelvinVec->getVariable();
        tmpVar1->setName( "TemperatureInKelvin" );
        meshAdapter->registerVectorAsData( TemperatureInKelvinVec );

        tmpVar1 = ResidualVec->getVariable();
        tmpVar1->setName( "Residual" );

        meshAdapter->registerVectorAsData( ResidualVec );
        manager->writeFile<AMP::Mesh::SiloIO>( exeName, 0 );
#endif
    }

    input_db.reset();

    ut.passes( exeName );

    AMP::AMPManager::shutdown();
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    node  = AMP::AMP_MPI::getRank();
    nodes = AMP::AMP_MPI::getNodes();

    try {
        linearThermalTest( ut );
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();
    } catch ( std::exception &err ) {
        std::cout << "ERROR: " << err.what() << std::endl;
        ut.numFails++;
    } catch ( ... ) {
        std::cout << "ERROR: "
                  << "An unknown exception was thrown." << std::endl;
        ut.numFails++;
    }
    return ut.numFails;
}
