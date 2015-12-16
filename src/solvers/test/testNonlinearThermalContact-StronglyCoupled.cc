#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include <string>

#include "materials/Material.h"
#include "utils/AMPManager.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Utilities.h"
#include "utils/shared_ptr.h"


#include "vectors/SimpleVector.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"

#include "ampmesh/MeshVariable.h"
#include "utils/Writer.h"

#include "operators/ColumnBoundaryOperator.h"
#include "operators/ElementOperationFactory.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/libmesh/MassLinearElement.h"
#include "operators/libmesh/MassLinearFEOperator.h"

#include "operators/CoupledOperator.h"
#include "operators/CoupledOperatorParameters.h"
#include "operators/DirichletMatrixCorrection.h"
#include "operators/DirichletVectorCorrection.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NeumannVectorCorrection.h"
#include "operators/NeumannVectorCorrectionParameters.h"
#include "operators/NeutronicsRhs.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/RobinMatrixCorrection.h"
#include "operators/RobinVectorCorrection.h"
#include "operators/libmesh/VolumeIntegralOperator.h"
#include "operators/map/MapOperatorParameters.h"
#include "operators/map/MapSurface.h"

#include "../ColumnSolver.h"
#include "../PetscKrylovSolver.h"
#include "../PetscKrylovSolverParameters.h"
#include "../PetscSNESSolver.h"
#include "../PetscSNESSolverParameters.h"
#include "../TrilinosMLSolver.h"


void thermalContactTest( AMP::UnitTest *ut, std::string exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    //  AMP::AMPManager::startup();
    //  AMP::Materials::initialize();

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP::PIO::logAllNodes( log_file );

    //  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
    //  std::string mesh_file = input_db->getString("Mesh");

    AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams(
        new AMP::Mesh::MeshManagerParameters( input_db ) );
    AMP::Mesh::MeshManager::shared_ptr manager( new AMP::Mesh::MeshManager( mgrParams ) );
    AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter1 = manager->getMesh( "pellet" );
    AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter2 = manager->getMesh( "clad" );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    AMP::LinearAlgebra::Variable::shared_ptr TemperatureVar(
        new AMP::Mesh::NodalScalarVariable( "Temperature" ) );
    AMP::LinearAlgebra::Variable::shared_ptr inputVariable1(
        new AMP::Mesh::NodalScalarVariable( "Temperature", meshAdapter1 ) );
    AMP::LinearAlgebra::Variable::shared_ptr inputVariable2(
        new AMP::Mesh::NodalScalarVariable( "Temperature", meshAdapter2 ) );

    AMP::LinearAlgebra::Variable::shared_ptr outputVariable1(
        new AMP::Mesh::NodalScalarVariable( "Temperature", meshAdapter1 ) );
    AMP::LinearAlgebra::Variable::shared_ptr outputVariable2(
        new AMP::Mesh::NodalScalarVariable( "Temperature", meshAdapter2 ) );

    double intguess = input_db->getDoubleWithDefault( "InitialGuess", 400 );

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvin =
        manager->createVector( TemperatureVar );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec =
        manager->createVector( TemperatureVar );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec = manager->createVector( TemperatureVar );

    TemperatureInKelvin->setToScalar( intguess );
    manager->registerVectorAsData( TemperatureInKelvin, "Temperature" );
    manager->registerVectorAsData( ResidualVec, "Residual" );


    //-----------------------------------------------
    //   CREATE THE NONLINEAR THERMAL OPERATOR 1 ----
    //-----------------------------------------------

    AMP_INSIST( input_db->keyExists( "NonlinearThermalOperator1" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel1;
    AMP::shared_ptr<AMP::Database> nonlinearThermalDatabase1 =
        input_db->getDatabase( "NonlinearThermalOperator1" );
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator1 =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter1, nonlinearThermalDatabase1, thermalTransportModel1 ) );

    // initialize the input variable
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator1 =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator1->getVolumeOperator() );

    // initialize the output variable
    // AMP::LinearAlgebra::Variable::shared_ptr outputVariable1 =
    // thermalVolumeOperator1->getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec1 =
        TemperatureInKelvin->subsetVectorForVariable( inputVariable1 );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec1 =
        RightHandSideVec->subsetVectorForVariable( outputVariable1 );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec1 =
        ResidualVec->subsetVectorForVariable( outputVariable1 );

    TemperatureInKelvinVec1->setToScalar( 400 );

    //-------------------------------------
    //   CREATE THE LINEAR THERMAL OPERATOR 1 ----
    //-------------------------------------

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel1;
    AMP::shared_ptr<AMP::InputDatabase> bvpDatabase1 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            input_db->getDatabase( "LinearThermalOperator1" ) );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator1 =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter1, bvpDatabase1, transportModel1 ) );

    //-------------------------------------
    //  CREATE THE NEUTRONICS SOURCE  //
    //-------------------------------------
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    AMP::shared_ptr<AMP::Database> neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    AMP::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(
        new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ) );
    neutronicsParams->d_MeshAdapter = meshAdapter1;
    AMP::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(
        new AMP::Operator::NeutronicsRhs( neutronicsParams ) );

    AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar =
        neutronicsOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr SpecificPowerVec =
        meshAdapter1->createVector( SpecificPowerVar );

    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    //----------------------------------------------------------
    //  Integrate Nuclear Rhs over Desnity * Volume //
    //----------------------------------------------------------

    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
    AMP::shared_ptr<AMP::Database> sourceDatabase =
        input_db->getDatabase( "VolumeIntegralOperator" );
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter1, sourceDatabase, stransportModel ) );

    // Create the power (heat source) vector.
    AMP::LinearAlgebra::Variable::shared_ptr PowerInWattsVar = sourceOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr PowerInWattsVec =
        meshAdapter1->createVector( PowerInWattsVar );
    PowerInWattsVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, PowerInWattsVec );

    // copy the power to pellet RHS vector
    RightHandSideVec1->copyVector( PowerInWattsVec );

    cout << "L2 Norm of the PowerInWattsVec and RightHSVec1" << PowerInWattsVec->L2Norm() << " "
         << RightHandSideVec1->L2Norm() << endl;
    //--------------------------------------

    AMP_INSIST( input_db->keyExists( "NonlinearSolver" ), "Key ''NonlinearSolver'' is missing!" );

    AMP::shared_ptr<AMP::Database> nonlinearSolver_db1 = input_db->getDatabase( "NonlinearSolver" );
    AMP::shared_ptr<AMP::Database> linearSolver_db1 =
        nonlinearSolver_db1->getDatabase( "LinearSolver" );

    //--------------------------------------------
    //   CREATE THE NONLINEAR THERMAL OPERATOR 2 ----
    //--------------------------------------------

    AMP_INSIST( input_db->keyExists( "NonlinearThermalOperator2" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel2;
    AMP::shared_ptr<AMP::Database> nonlinearThermalDatabase2 =
        input_db->getDatabase( "NonlinearThermalOperator2" );
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator2 =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter2, nonlinearThermalDatabase2, thermalTransportModel2 ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator2 =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator2->getVolumeOperator() );

    // initialize the output variable
    // AMP::LinearAlgebra::Variable::shared_ptr outputVariable2 =
    // thermalVolumeOperator2->getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec2 =
        TemperatureInKelvin->subsetVectorForVariable( inputVariable2 );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec2 =
        RightHandSideVec->subsetVectorForVariable( outputVariable2 );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec2 =
        ResidualVec->subsetVectorForVariable( outputVariable2 );

    TemperatureInKelvinVec2->setToScalar( 400 );
    //--------------------------------------------
    //   CREATE THE LINEAR THERMAL OPERATOR 2 ----
    //--------------------------------------------

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel2;
    AMP::shared_ptr<AMP::InputDatabase> bvpDatabase2 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            input_db->getDatabase( "LinearThermalOperator2" ) );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator2 =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter2, bvpDatabase2, transportModel2 ) );

    //-------------------------------------
    AMP::shared_ptr<AMP::InputDatabase> mapcladtopellet_db =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase( "MapCladtoPellet" ) );
    AMP::shared_ptr<AMP::Operator::MapSurface> mapcladtopellet =
        AMP::dynamic_pointer_cast<AMP::Operator::MapSurface>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter2, meshAdapter1, mapcladtopellet_db ) );

    AMP::shared_ptr<AMP::InputDatabase> mappellettoclad_db =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase( "MapPellettoClad" ) );
    AMP::shared_ptr<AMP::Operator::MapSurface> mappellettoclad =
        AMP::dynamic_pointer_cast<AMP::Operator::MapSurface>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter1, meshAdapter2, mappellettoclad_db ) );

    //------------------------------------------

    AMP::Operator::Operator::shared_ptr boundaryOp1;
    boundaryOp1 = nonlinearThermalOperator1->getBoundaryOperator();

    AMP::Operator::Operator::shared_ptr robinBoundaryOp1;
    robinBoundaryOp1 =
        ( AMP::dynamic_pointer_cast<AMP::Operator::BoundaryOperator>( boundaryOp1 ) );

    AMP::shared_ptr<AMP::InputDatabase> boundaryDatabase1 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            nonlinearThermalDatabase1->getDatabase( "BoundaryOperator" ) );
    AMP::shared_ptr<AMP::InputDatabase> robinboundaryDatabase1 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( boundaryDatabase1 );

    robinboundaryDatabase1->putBool( "constant_flux", false );
    robinboundaryDatabase1->putBool( "skip_matrix_correction", true );
    AMP::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters1(
        new AMP::Operator::NeumannVectorCorrectionParameters( robinboundaryDatabase1 ) );


    //------------------------------------------

    AMP::Operator::Operator::shared_ptr boundaryOp2;
    boundaryOp2 = nonlinearThermalOperator2->getBoundaryOperator();

    AMP::Operator::Operator::shared_ptr robinBoundaryOp2;
    robinBoundaryOp2 =
        ( AMP::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator>( boundaryOp2 ) )
            ->getBoundaryOperator( 0 );

    AMP::shared_ptr<AMP::InputDatabase> boundaryDatabase2 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            nonlinearThermalDatabase2->getDatabase( "BoundaryOperator" ) );
    AMP::shared_ptr<AMP::InputDatabase> robinboundaryDatabase2 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            boundaryDatabase2->getDatabase( "RobinVectorCorrection" ) );

    robinboundaryDatabase2->putBool( "constant_flux", false );
    robinboundaryDatabase2->putBool( "skip_matrix_correction", true );
    AMP::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters2(
        new AMP::Operator::NeumannVectorCorrectionParameters( robinboundaryDatabase2 ) );


    //-------------------------------------
    // Coupling Map to the Nonlinear Operators
    AMP::shared_ptr<AMP::InputDatabase> tmp_db( new AMP::InputDatabase( "Dummy" ) );
    AMP::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledNonlinearPelletParams(
        new AMP::Operator::CoupledOperatorParameters( tmp_db ) );

    coupledNonlinearPelletParams->d_MapOperator = mapcladtopellet;
    coupledNonlinearPelletParams->d_BVPOperator = nonlinearThermalOperator1;
    AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledNonlinearPellet(
        new AMP::Operator::CoupledOperator( coupledNonlinearPelletParams ) );
    //-------------------------------------
    AMP::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledNonlinearCladParams(
        new AMP::Operator::CoupledOperatorParameters( tmp_db ) );
    coupledNonlinearCladParams->d_MapOperator = mappellettoclad;
    coupledNonlinearCladParams->d_BVPOperator = nonlinearThermalOperator2;
    AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledNonlinearClad(
        new AMP::Operator::CoupledOperator( coupledNonlinearCladParams ) );

    //-------------------------------------
    // Column of Coupled Operators
    AMP::shared_ptr<AMP::Operator::OperatorParameters> nonlinearParams(
        new AMP::Operator::OperatorParameters( tmp_db ) );
    AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearCoupledOperator(
        new AMP::Operator::ColumnOperator( nonlinearParams ) );
    nonlinearCoupledOperator->append( coupledNonlinearPellet );
    nonlinearCoupledOperator->append( coupledNonlinearClad );

    /*-------------------------------------
      //Coupling Map to the Linear Operators
      AMP::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledLinearPelletParams(new
      AMP::Operator::CoupledOperatorParameters(tmp_db));
      coupledLinearPelletParams->d_MapOperator = mapcladtopellet ;
      coupledLinearPelletParams->d_BVPOperator = linearThermalOperator1;
      AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledLinearPellet(new
      AMP::Operator::CoupledOperator(coupledLinearPelletParams));
      //-------------------------------------
      AMP::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledLinearCladParams(new
      AMP::Operator::CoupledOperatorParameters(tmp_db));
      coupledLinearCladParams->d_MapOperator = mappellettoclad ;
      coupledLinearCladParams->d_BVPOperator = linearThermalOperator2;
      AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledLinearClad(new
      AMP::Operator::CoupledOperator(coupledLinearCladParams));
    */
    //-------------------------------------
    // Column of Coupled Operators
    AMP::shared_ptr<AMP::Operator::OperatorParameters> linearParams(
        new AMP::Operator::OperatorParameters( tmp_db ) );
    AMP::shared_ptr<AMP::Operator::ColumnOperator> linearCoupledOperator(
        new AMP::Operator::ColumnOperator( linearParams ) );
    linearCoupledOperator->append( linearThermalOperator1 );
    linearCoupledOperator->append( linearThermalOperator2 );

    //------------------------------------------------------------------

    AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    AMP::shared_ptr<AMP::Database> linearSolver_db =
        nonlinearSolver_db->getDatabase( "LinearSolver" );

    //----------------------------------------------------------------//
    // initialize the nonlinear solver
    AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(
        new AMP::Solver::PetscSNESSolverParameters( nonlinearSolver_db ) );

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = nonlinearCoupledOperator;
    nonlinearSolverParams->d_pInitialGuess = TemperatureInKelvin;
    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(
        new AMP::Solver::PetscSNESSolver( nonlinearSolverParams ) );

    //-------------------------------------------------------------------------//
    // initialize the column preconditioner which is a diagonal block preconditioner
    AMP::shared_ptr<AMP::Database> columnPreconditioner_db =
        linearSolver_db->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> columnPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( columnPreconditioner_db ) );
    columnPreconditionerParams->d_pOperator = linearCoupledOperator;
    AMP::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(
        new AMP::Solver::ColumnSolver( columnPreconditionerParams ) );

    AMP::shared_ptr<AMP::Database> pelletPreconditioner_db =
        columnPreconditioner_db->getDatabase( "pelletPreconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> pelletPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( pelletPreconditioner_db ) );
    pelletPreconditionerParams->d_pOperator = linearThermalOperator1;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearPelletPreconditioner(
        new AMP::Solver::TrilinosMLSolver( pelletPreconditionerParams ) );

    AMP::shared_ptr<AMP::Database> cladPreconditioner_db =
        columnPreconditioner_db->getDatabase( "cladPreconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> cladPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( cladPreconditioner_db ) );
    cladPreconditionerParams->d_pOperator = linearThermalOperator2;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearCladPreconditioner(
        new AMP::Solver::TrilinosMLSolver( cladPreconditionerParams ) );

    columnPreconditioner->append( linearPelletPreconditioner );
    columnPreconditioner->append( linearCladPreconditioner );

    //--------------------------------------------------------------------//
    // register the preconditioner with the Jacobian free Krylov solver
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver =
        nonlinearSolver->getKrylovSolver();

    linearSolver->setPreconditioner( columnPreconditioner );

    //--------------------------------------------------------------------//

    AMP::LinearAlgebra::Variable::shared_ptr TemperatureMapVar(
        new AMP::Mesh::NodalScalarVariable( "TemperatureMap" ) );
    AMP::LinearAlgebra::Variable::shared_ptr inputMapVariable1(
        new AMP::Mesh::NodalScalarVariable( "TemperatureMap", meshAdapter1 ) );
    AMP::LinearAlgebra::Variable::shared_ptr inputMapVariable2(
        new AMP::Mesh::NodalScalarVariable( "TemperatureMap", meshAdapter2 ) );

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureMapvector =
        manager->createVector( TemperatureMapVar );
    AMP::LinearAlgebra::Vector::shared_ptr scratchCladToPellet =
        TemperatureMapvector->subsetVectorForVariable( inputMapVariable1 );
    AMP::LinearAlgebra::Vector::shared_ptr scratchPelletToClad =
        TemperatureMapvector->subsetVectorForVariable( inputMapVariable2 );

    correctionParameters1->d_variableFlux = scratchCladToPellet;
    correctionParameters2->d_variableFlux = scratchPelletToClad;

    robinBoundaryOp1->reset( correctionParameters1 );
    robinBoundaryOp2->reset( correctionParameters2 );

    mapcladtopellet->setVector( scratchCladToPellet );
    mappellettoclad->setVector( scratchPelletToClad );
    //-------------------------------------
    // Applying boundary conditions to the nonlinear BVP Operator

    nonlinearThermalOperator1->modifyRHSvector( RightHandSideVec1 );
    nonlinearThermalOperator1->modifyInitialSolutionVector( TemperatureInKelvinVec1 );

    nonlinearThermalOperator2->modifyRHSvector( RightHandSideVec2 );
    nonlinearThermalOperator2->modifyInitialSolutionVector( TemperatureInKelvinVec2 );
    //-------------------------------------

    nonlinearCoupledOperator->residual( RightHandSideVec, TemperatureInKelvin, ResidualVec );

    double initialResidualNorm = ResidualVec->L2Norm();

    AMP::pout << "Initial Residual Norm: " << initialResidualNorm << std::endl;

    nonlinearSolver->setZeroInitialGuess( false );

    nonlinearSolver->solve( RightHandSideVec, TemperatureInKelvin );

    bool testPassed = false;

    double finalResidualNorm = ResidualVec->L2Norm();

    AMP::pout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    if ( finalResidualNorm < 1.e-5 ) {
        testPassed = true;
    }

    nonlinearCoupledOperator->residual( RightHandSideVec, TemperatureInKelvin, ResidualVec );

#ifdef USE_EXT_SILO
    manager->writeFile<AMP::SiloIO>( exeName, 0 );
#endif

    //-------------------------------------

    if ( testPassed ) {
        ut.passes( "Coupled solve of Composite Operator using simple Map and Robin Operators. " );
    } else {
        ITFAILS;
    }

    input_db.reset();

    ut.passes( exeName );

    //  AMP::AMPManager::shutdown();
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    try {
        thermalContactTest( ut, "testNonlinearThermalContactCoupled_HALDEN" );
    } catch ( std::exception &err ) {
        std::cout << "ERROR: While testing " << argv[0] << err.what() << std::endl;
        ut.failure( "ERROR: While testing" );
    } catch ( ... ) {
        std::cout << "ERROR: While testing " << argv[0] << "An unknown exception was thrown."
                  << std::endl;
        ut.failure( "ERROR: While testing" );
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
