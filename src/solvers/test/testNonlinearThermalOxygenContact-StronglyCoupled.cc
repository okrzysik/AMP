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
    AMP::LinearAlgebra::Variable::shared_ptr inputThermalVariable1(
        new AMP::Mesh::NodalScalarVariable( "Temperature", meshAdapter1 ) );
    AMP::LinearAlgebra::Variable::shared_ptr inputThermalVariable2(
        new AMP::Mesh::NodalScalarVariable( "Temperature", meshAdapter2 ) );

    AMP::LinearAlgebra::Variable::shared_ptr outputThermalVariable1(
        new AMP::Mesh::NodalScalarVariable( "Temperature", meshAdapter1 ) );
    AMP::LinearAlgebra::Variable::shared_ptr outputThermalVariable2(
        new AMP::Mesh::NodalScalarVariable( "Temperature", meshAdapter2 ) );

    AMP::LinearAlgebra::Variable::shared_ptr ConcentrationVar(
        new AMP::Mesh::NodalScalarVariable( "Concentration" ) );
    AMP::LinearAlgebra::Variable::shared_ptr inputOxygenVariable1(
        new AMP::Mesh::NodalScalarVariable( "Concentration", meshAdapter1 ) );
    AMP::LinearAlgebra::Variable::shared_ptr outputOxygenVariable1(
        new AMP::Mesh::NodalScalarVariable( "Concentration", meshAdapter1 ) );

    AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> TempOxyVariable(
        new AMP::LinearAlgebra::MultiVariable( "inputVariable" ) );
    TempOxyVariable->add( inputThermalVariable1 );
    TempOxyVariable->add( inputThermalVariable2 );
    TempOxyVariable->add( inputOxygenVariable1 );

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureOxygenSolution =
        manager->createVector( TempOxyVariable );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec =
        manager->createVector( TempOxyVariable );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec = manager->createVector( TempOxyVariable );

    //  manager->registerVectorAsData ( TemperatureOxygenSolution , "Solution" );
    //  manager->registerVectorAsData ( ResidualVec               , "Residual" );

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

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureSolutionVec1 =
        TemperatureOxygenSolution->subsetVectorForVariable( inputThermalVariable1 );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec1 =
        RightHandSideVec->subsetVectorForVariable( outputThermalVariable1 );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec1 =
        ResidualVec->subsetVectorForVariable( outputThermalVariable1 );

    TemperatureSolutionVec1->setToScalar( 400 );

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

    neutronicsOperator->apply( nullVec, nullVec, SpecificPowerVec, 1., 0. );

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
    sourceOperator->apply( nullVec, SpecificPowerVec, PowerInWattsVec, 1., 0. );

    // copy the power to pellet RHS vector
    RightHandSideVec1->copyVector( PowerInWattsVec );

    cout << "L2 Norm of the PowerInWattsVec and RightHSVec1" << PowerInWattsVec->L2Norm() << " "
         << RightHandSideVec1->L2Norm() << endl;

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear oxygen diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearOxygenOperator1" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> oxygenTransportModel;
    AMP::shared_ptr<AMP::Database> nonlinearOxygenDatabase =
        input_db->getDatabase( "testNonlinearOxygenOperator1" );
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearOxygenOperator1 =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter1, nonlinearOxygenDatabase, oxygenTransportModel ) );

    AMP::LinearAlgebra::Vector::shared_ptr ConcentrationSolutionVec1 =
        TemperatureOxygenSolution->subsetVectorForVariable( inputOxygenVariable1 );
    AMP::LinearAlgebra::Vector::shared_ptr ConcentrationRightHandSideVec1 =
        RightHandSideVec->subsetVectorForVariable( outputOxygenVariable1 );
    AMP::LinearAlgebra::Vector::shared_ptr ConcentrationResidualVec1 =
        ResidualVec->subsetVectorForVariable( outputOxygenVariable1 );

    ConcentrationSolutionVec1->setToScalar( 0.01 );
    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for oxygen
    AMP_INSIST( input_db->keyExists( "testLinearOxygenOperator1" ), "key missing!" );
    AMP::shared_ptr<AMP::Database> linearOxygenDatabase =
        input_db->getDatabase( "testLinearOxygenOperator1" );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearOxygenOperator1 =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter1, linearOxygenDatabase, oxygenTransportModel ) );

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
    // AMP::LinearAlgebra::Variable::shared_ptr outputThermalVariable2 =
    // thermalVolumeOperator2->getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureSolutionVec2 =
        TemperatureOxygenSolution->subsetVectorForVariable( inputThermalVariable2 );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec2 =
        RightHandSideVec->subsetVectorForVariable( outputThermalVariable2 );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec2 =
        ResidualVec->subsetVectorForVariable( outputThermalVariable2 );

    TemperatureSolutionVec2->setToScalar( 400 );
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

    //--------------------------------------

    AMP_INSIST( input_db->keyExists( "NonlinearSolver" ), "Key ''NonlinearSolver'' is missing!" );

    AMP::shared_ptr<AMP::Database> nonlinearSolver_db1 = input_db->getDatabase( "NonlinearSolver" );
    AMP::shared_ptr<AMP::Database> linearSolver_db1 =
        nonlinearSolver_db1->getDatabase( "LinearSolver" );


    //-------------------------------------
    // Coupling Map to the Nonlinear Operators
    /*
      // create a column operator object for nonlinear thermal-oxygen diffusion on
      // Pellet
      AMP::shared_ptr<AMP::Operator::OperatorParameters> params;
      AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearThermalOxygenOperator1(new
      AMP::Operator::ColumnOperator(params));
      nonlinearThermalOxygenOperator1->append(nonlinearOxygenOperator1);
      nonlinearThermalOxygenOperator1->append(nonlinearThermalOperator1);
    */
    //-------------------------------------

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
    nonlinearCoupledOperator->append( nonlinearOxygenOperator1 );

    //---------------------------------------------------------------
    /*
      // create a column operator object for linear thermal-oxygen
      AMP::shared_ptr<AMP::Operator::ColumnOperator> linearThermalOxygenOperator1(new
      AMP::Operator::ColumnOperator(params));
      linearThermalOxygenOperator1->append(linearOxygenOperator1);
      linearThermalOxygenOperator1->append(linearThermalOperator1);
    */
    //---------------------------------------
    // Column of Coupled Operators
    AMP::shared_ptr<AMP::Operator::OperatorParameters> linearParams(
        new AMP::Operator::OperatorParameters( tmp_db ) );
    AMP::shared_ptr<AMP::Operator::ColumnOperator> linearCoupledOperator(
        new AMP::Operator::ColumnOperator( linearParams ) );
    linearCoupledOperator->append( linearThermalOperator1 );
    linearCoupledOperator->append( linearThermalOperator2 );
    linearCoupledOperator->append( linearOxygenOperator1 );

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
    nonlinearSolverParams->d_pInitialGuess = TemperatureOxygenSolution;
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

    AMP::shared_ptr<AMP::Database> pelletThermalPreconditioner_db =
        columnPreconditioner_db->getDatabase( "pelletThermalPreconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> pelletThermalPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( pelletThermalPreconditioner_db ) );
    pelletThermalPreconditionerParams->d_pOperator = linearThermalOperator1;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearPelletThermalPreconditioner(
        new AMP::Solver::TrilinosMLSolver( pelletThermalPreconditionerParams ) );

    AMP::shared_ptr<AMP::Database> pelletOxygenPreconditioner_db =
        columnPreconditioner_db->getDatabase( "pelletOxygenPreconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> pelletOxygenPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( pelletOxygenPreconditioner_db ) );
    pelletOxygenPreconditionerParams->d_pOperator = linearOxygenOperator1;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearPelletOxygenPreconditioner(
        new AMP::Solver::TrilinosMLSolver( pelletOxygenPreconditionerParams ) );

    AMP::shared_ptr<AMP::Database> cladThermalPreconditioner_db =
        columnPreconditioner_db->getDatabase( "cladThermalPreconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> cladThermalPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( cladThermalPreconditioner_db ) );
    cladThermalPreconditionerParams->d_pOperator = linearThermalOperator2;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearCladThermalPreconditioner(
        new AMP::Solver::TrilinosMLSolver( cladThermalPreconditionerParams ) );

    columnPreconditioner->append( linearPelletThermalPreconditioner );
    columnPreconditioner->append( linearCladThermalPreconditioner );
    columnPreconditioner->append( linearPelletOxygenPreconditioner );

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
    nonlinearThermalOperator1->modifyInitialSolutionVector( TemperatureSolutionVec1 );

    nonlinearThermalOperator2->modifyRHSvector( RightHandSideVec2 );
    nonlinearThermalOperator2->modifyInitialSolutionVector( TemperatureSolutionVec2 );

    nonlinearOxygenOperator1->modifyRHSvector( ConcentrationRightHandSideVec1 );
    nonlinearOxygenOperator1->modifyInitialSolutionVector( ConcentrationSolutionVec1 );

    //-------------------------------------
    /*
      nonlinearThermalOxygenOperator1->apply(nullVec, TemperatureOxygenSolution, ResidualVec, 1.0,
      0.0);
      linearThermalOxygenOperator1->reset(nonlinearThermalOxygenOperator1->getParameters("Jacobian",
      TemperatureOxygenSolution));
    */
    //-------------------------------------

    nonlinearCoupledOperator->apply(
        RightHandSideVec, TemperatureOxygenSolution, ResidualVec, 1.0, -1.0 );

    double initialResidualNorm = ResidualVec->L2Norm();

    AMP::pout << "Initial Residual Norm: " << initialResidualNorm << std::endl;

    nonlinearSolver->setZeroInitialGuess( false );

    nonlinearSolver->solve( RightHandSideVec, TemperatureOxygenSolution );

    bool testPassed = false;

    double finalResidualNorm = ResidualVec->L2Norm();

    AMP::pout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    if ( finalResidualNorm < 1.e-5 ) {
        testPassed = true;
    }

    nonlinearCoupledOperator->apply(
        RightHandSideVec, TemperatureOxygenSolution, ResidualVec, 1.0, -1.0 );

    meshAdapter1->registerVectorAsData( TemperatureSolutionVec1, "Temperature_pellet" );
    meshAdapter2->registerVectorAsData( TemperatureSolutionVec2, "Temperature_clad" );
    meshAdapter1->registerVectorAsData( ConcentrationSolutionVec1, "Concentration_pellet" );

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

    thermalContactTest( ut, "testNonlinearThermalOxygenContactCoupled_HALDEN" );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
