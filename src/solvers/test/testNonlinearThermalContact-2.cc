#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <string>

#include "materials/Material.h"
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
#include "operators/map/Map1Dto3D.h"
#include "operators/map/Map3Dto1D.h"
#include "operators/map/MapOperatorParameters.h"

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

    double intguess = input_db->getDoubleWithDefault( "InitialGuess", 400 );

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvin =
        manager->createVector( TemperatureVar );
    TemperatureInKelvin->setToScalar( intguess );
    manager->registerVectorAsData( TemperatureInKelvin, "Temperature" );


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
    AMP::LinearAlgebra::Variable::shared_ptr outputVariable1 =
        thermalVolumeOperator1->getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec1 =
        TemperatureInKelvin->subsetVectorForVariable( inputVariable1 );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec1 =
        meshAdapter1->createVector( outputVariable1 );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec1 =
        meshAdapter1->createVector( outputVariable1 );

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
    //  Integrate Nuclear Rhs over Desnity * GeomType::Volume //
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

    //--------------------------------------

    AMP_INSIST( input_db->keyExists( "NonlinearSolver" ), "Key ''NonlinearSolver'' is missing!" );

    AMP::shared_ptr<AMP::Database> nonlinearSolver_db1 = input_db->getDatabase( "NonlinearSolver" );
    AMP::shared_ptr<AMP::Database> linearSolver_db1 =
        nonlinearSolver_db1->getDatabase( "LinearSolver" );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the nonlinear solver
    AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams1(
        new AMP::Solver::PetscSNESSolverParameters( nonlinearSolver_db1 ) );

    // change the next line to get the correct communicator out
    nonlinearSolverParams1->d_comm          = globalComm;
    nonlinearSolverParams1->d_pOperator     = nonlinearThermalOperator1;
    nonlinearSolverParams1->d_pInitialGuess = TemperatureInKelvinVec1;

    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver1(
        new AMP::Solver::PetscSNESSolver( nonlinearSolverParams1 ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//

    AMP::shared_ptr<AMP::Database> thermalPreconditioner_db1 =
        linearSolver_db1->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams1(
        new AMP::Solver::SolverStrategyParameters( thermalPreconditioner_db1 ) );
    thermalPreconditionerParams1->d_pOperator = linearThermalOperator1;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearThermalPreconditioner1(
        new AMP::Solver::TrilinosMLSolver( thermalPreconditionerParams1 ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // register the preconditioner with the Jacobian free Krylov solver
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver1 =
        nonlinearSolver1->getKrylovSolver();
    linearSolver1->setPreconditioner( linearThermalPreconditioner1 );
    nonlinearThermalOperator1->residual( RightHandSideVec1, TemperatureInKelvinVec1, ResidualVec1 );

    //---------------------------------------------
    //     CREATE THE CONTACT GAP OPERATOR
    //---------------------------------------------

    AMP_INSIST( input_db->keyExists( "GapOperator" ), "Key ''GapOperator'' is missing!" );
    AMP::shared_ptr<AMP::InputDatabase> gapDatabase =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase( "GapOperator" ) );

    double heff = ( gapDatabase )->getDouble( "Convective_Coefficient" );
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> gapVariable(
        new AMP::LinearAlgebra::Variable( "Gap" ) );


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

    // initialize the input variable
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator2 =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator2->getVolumeOperator() );

    // initialize the output variable
    AMP::LinearAlgebra::Variable::shared_ptr outputVariable2 =
        thermalVolumeOperator2->getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec2 =
        TemperatureInKelvin->subsetVectorForVariable( inputVariable2 );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec2 =
        meshAdapter2->createVector( outputVariable2 );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec2 =
        meshAdapter2->createVector( outputVariable2 );

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

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the nonlinear solver
    AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams2(
        new AMP::Solver::PetscSNESSolverParameters( nonlinearSolver_db1 ) );

    // change the next line to get the correct communicator out
    nonlinearSolverParams2->d_comm          = globalComm;
    nonlinearSolverParams2->d_pOperator     = nonlinearThermalOperator2;
    nonlinearSolverParams2->d_pInitialGuess = TemperatureInKelvinVec2;

    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver2(
        new AMP::Solver::PetscSNESSolver( nonlinearSolverParams2 ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//

    AMP::shared_ptr<AMP::Database> thermalPreconditioner_db2 =
        linearSolver_db1->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams2(
        new AMP::Solver::SolverStrategyParameters( thermalPreconditioner_db1 ) );
    thermalPreconditionerParams2->d_pOperator = linearThermalOperator2;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearThermalPreconditioner2(
        new AMP::Solver::TrilinosMLSolver( thermalPreconditionerParams2 ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // register the preconditioner with the Jacobian free Krylov solver
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver2 =
        nonlinearSolver2->getKrylovSolver();
    linearSolver2->setPreconditioner( linearThermalPreconditioner2 );
    nonlinearThermalOperator2->residual( RightHandSideVec2, TemperatureInKelvinVec2, ResidualVec2 );

    //-------------------------------------

    AMP::LinearAlgebra::Vector::shared_ptr variableFluxVec1 =
        meshAdapter1->createVector( inputVariable1 );
    AMP::LinearAlgebra::Vector::shared_ptr scratchTempVec1 =
        meshAdapter1->createVector( inputVariable1 );
    variableFluxVec1->setToScalar( 0.0 );

    AMP::LinearAlgebra::Vector::shared_ptr variableFluxVec2 =
        meshAdapter2->createVector( inputVariable2 );
    AMP::LinearAlgebra::Vector::shared_ptr scratchTempVec2 =
        meshAdapter2->createVector( inputVariable2 );
    variableFluxVec2->setToScalar( 0.0 );

    //-------------------------------------

    AMP::shared_ptr<AMP::InputDatabase> map3dto1d_db1 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase( "MapPelletto1D" ) );
    AMP::shared_ptr<AMP::Operator::MapOperatorParameters> map3dto1dParams1(
        new AMP::Operator::MapOperatorParameters( map3dto1d_db1 ) );
    map3dto1dParams1->d_MeshAdapter = meshAdapter1;
    AMP::shared_ptr<AMP::Operator::Map3Dto1D> map1ToLowDim(
        new AMP::Operator::Map3Dto1D( map3dto1dParams1 ) );

    AMP::shared_ptr<AMP::InputDatabase> map1dto3d_db1 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase( "Map1DtoClad" ) );
    AMP::shared_ptr<AMP::Operator::MapOperatorParameters> map1dto3dParams1(
        new AMP::Operator::MapOperatorParameters( map1dto3d_db1 ) );
    map1dto3dParams1->d_MapAdapter = meshAdapter2;
    //-------------------------------------
    // This is related to But # 1219 and 1210.
    //  -- It dies in compute_Z_locations of the constructor for mat1dto3d.
    ut.passes( "Everything up till constructing 1Dto3D passes." );
    // if( AMP::AMP_MPI::getNodes() == 1 ) {
    //-------------------------------------
    AMP::shared_ptr<AMP::Operator::Map1Dto3D> map1ToHighDim(
        new AMP::Operator::Map1Dto3D( map1dto3dParams1 ) );

    map1ToLowDim->setZLocations( map1ToHighDim->getZLocations() );

    AMP::shared_ptr<AMP::InputDatabase> map3dto1d_db2 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase( "MapCladto1D" ) );
    AMP::shared_ptr<AMP::Operator::MapOperatorParameters> map3dto1dParams2(
        new AMP::Operator::MapOperatorParameters( map3dto1d_db2 ) );
    map3dto1dParams2->d_MeshAdapter = meshAdapter2;
    AMP::shared_ptr<AMP::Operator::Map3Dto1D> map2ToLowDim(
        new AMP::Operator::Map3Dto1D( map3dto1dParams2 ) );

    AMP::shared_ptr<AMP::InputDatabase> map1dto3d_db2 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase( "Map1DtoPellet" ) );
    AMP::shared_ptr<AMP::Operator::MapOperatorParameters> map1dto3dParams2(
        new AMP::Operator::MapOperatorParameters( map1dto3d_db2 ) );
    map1dto3dParams2->d_MapAdapter = meshAdapter1;
    AMP::shared_ptr<AMP::Operator::Map1Dto3D> map2ToHighDim(
        new AMP::Operator::Map1Dto3D( map1dto3dParams2 ) );

    map2ToLowDim->setZLocations( map2ToHighDim->getZLocations() );

    //------------------------------------------

    AMP::Operator::Operator::shared_ptr boundaryOp1;
    boundaryOp1 = nonlinearThermalOperator1->getBoundaryOperator();

    AMP::Operator::Operator::shared_ptr robinBoundaryOp1;
    //  robinBoundaryOp1 =
    //  (AMP::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator>(boundaryOp1)
    //  )->getBoundaryOperator(0);
    robinBoundaryOp1 =
        ( AMP::dynamic_pointer_cast<AMP::Operator::BoundaryOperator>( boundaryOp1 ) );

    AMP::shared_ptr<AMP::InputDatabase> boundaryDatabase1 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            nonlinearThermalDatabase1->getDatabase( "BoundaryOperator" ) );
    //  AMP::shared_ptr<AMP::InputDatabase> robinboundaryDatabase1 =
    //  AMP::dynamic_pointer_cast<AMP::InputDatabase>(
    //  boundaryDatabase1->getDatabase("RobinVectorCorrection"));
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

    size_t gapVecCladSize = map1ToHighDim->getNumZlocations();
    AMP::LinearAlgebra::Vector::shared_ptr gapVecClad =
        AMP::LinearAlgebra::SimpleVector<double>::create( gapVecCladSize, gapVariable );

    size_t gapVecPelletSize = map2ToHighDim->getNumZlocations();
    AMP::LinearAlgebra::Vector::shared_ptr gapVecPellet =
        AMP::LinearAlgebra::SimpleVector<double>::create( gapVecPelletSize, gapVariable );

    int cnt                                        = 0;
    AMP::LinearAlgebra::Vector::shared_ptr vecLag1 = meshAdapter1->createVector( outputVariable1 );
    vecLag1->copyVector( TemperatureInKelvinVec1 );
    AMP::LinearAlgebra::Vector::shared_ptr vecLag2 = meshAdapter2->createVector( outputVariable2 );
    vecLag2->copyVector( TemperatureInKelvinVec2 );

    bool testPassed = false;

    int maxIt = input_db->getIntegerWithDefault( "max_iterations", 100 );

    while ( cnt < maxIt ) {
        cnt++;

        RightHandSideVec1->zero();
        RightHandSideVec2->zero();

        RightHandSideVec1->copyVector( PowerInWattsVec );
        std::cout << "PowerInWattsVec norm  inside loop = " << RightHandSideVec1->L2Norm() << "\n";

        map2ToLowDim->apply( TemperatureInKelvinVec2, gapVecPellet );
        map2ToHighDim->apply( gapVecPellet, scratchTempVec1 );

        scratchTempVec1->scale( heff );
        variableFluxVec1->copyVector( scratchTempVec1 );

        correctionParameters1->d_variableFlux = variableFluxVec1;
        robinBoundaryOp1->reset( correctionParameters1 );

        std::cout << "Variable flux1 norm inside loop : " << variableFluxVec1->L2Norm() << endl;

        nonlinearThermalOperator1->modifyRHSvector( RightHandSideVec1 );
        nonlinearThermalOperator1->modifyInitialSolutionVector( TemperatureInKelvinVec1 );
        nonlinearSolver1->solve( RightHandSideVec1, TemperatureInKelvinVec1 );
        nonlinearThermalOperator1->residual(
            RightHandSideVec1, TemperatureInKelvinVec1, ResidualVec1 );

        std::cout << "Norm of TemperatureInKelvinVec1: " << TemperatureInKelvinVec1->L2Norm()
                  << endl;

        //------------------------------------------------------------
        map1ToLowDim->apply( TemperatureInKelvinVec1, gapVecClad );
        map1ToHighDim->apply( gapVecClad, scratchTempVec2 );

        scratchTempVec2->scale( heff );
        variableFluxVec2->copyVector( scratchTempVec2 );

        correctionParameters2->d_variableFlux = variableFluxVec2;
        robinBoundaryOp2->reset( correctionParameters2 );

        std::cout << "Variable flux2 norm inside loop : " << variableFluxVec2->L2Norm() << endl;

        nonlinearThermalOperator2->modifyRHSvector( RightHandSideVec2 );
        nonlinearThermalOperator2->modifyInitialSolutionVector( TemperatureInKelvinVec2 );
        nonlinearThermalOperator2->residual(
            RightHandSideVec2, TemperatureInKelvinVec2, ResidualVec2 );
        nonlinearSolver2->solve( RightHandSideVec2, TemperatureInKelvinVec2 );

        std::cout << "Residual Norm on Pellet after " << cnt
                  << " iteration is : " << ResidualVec1->L2Norm() << std::endl;
        std::cout << "Residual Norm on Clad after " << cnt
                  << " iteration is : " << ResidualVec2->L2Norm() << std::endl;

        vecLag2->subtract( TemperatureInKelvinVec2, vecLag2 );

//          if( nodes == 2 ) {
#ifdef USE_EXT_SILO
        manager->writeFile<AMP::Mesh::SiloIO>( exeName, 0 );
#endif
        //          }
        if ( vecLag2->L2Norm() < 1.e-6 ) {
            testPassed = true;
            break;
        } else {
            std::cout << "for iteration cnt = " << cnt << " --> " << vecLag1->L2Norm() << " "
                      << vecLag2->L2Norm() << std::endl;
        }
        std::cout << std::endl;

        vecLag1->copyVector( TemperatureInKelvinVec1 );
        vecLag2->copyVector( TemperatureInKelvinVec2 );
    }

    //-------------------------------------

    if ( testPassed ) {
        ut.passes( "Seggregated solve of Composite Operator using control loop of Nonlinear "
                   "Thermal+Robin->Map->Gap->Map->Ninlinear Thermal+Robin ." );
    } else {
        ITFAILS;
    }

    //} else {
    //  ut.expected_failure("parallel map3D-1D and map1D-3D fail in parallel, see bug #1219.");
    //}
    input_db.reset();

    ut.passes( exeName );

    //  AMP::AMPManager::shutdown();
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    thermalContactTest( ut, "testNonlinearThermalContactPicard2_HALDEN" );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
