#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/LinearOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"

#include "../applyTests.h"

#include <iostream>
#include <string>


static void thermoMechanicsTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db    = input_db->getDatabase( "Mesh" );
    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    // create a nonlinear BVP operator for nonlinear mechanics
    AMP_INSIST( input_db->keyExists( "testNonlinearMechanicsOperator" ), "key missing!" );

    auto nonlinearMechanicsOperator =
        std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearMechanicsOperator", input_db ) );

    auto nonlinearMechanicsVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsOperator->getVolumeOperator() );
    auto mechanicsMaterialModel = nonlinearMechanicsVolumeOperator->getMaterialModel();

    // create a nonlinear BVP operator for nonlinear thermal diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearThermalOperator" ), "key missing!" );

    auto nonlinearThermalOperator = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testNonlinearThermalOperator", input_db ) );

    auto nonlinearThermalVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator->getVolumeOperator() );
    auto thermalTransportModel = nonlinearThermalVolumeOperator->getTransportModel();

    // create a column operator object for nonlinear thermomechanics
    auto nonlinearThermoMechanicsOperator = std::make_shared<AMP::Operator::ColumnOperator>();
    nonlinearThermoMechanicsOperator->append( nonlinearMechanicsOperator );
    nonlinearThermoMechanicsOperator->append( nonlinearThermalOperator );


    // Create the relavent DOF managers
    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    auto nodalDofMap    = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
    int displacementDOFsPerNode = 3;
    auto displDofMap            = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, displacementDOFsPerNode, split );

    // initialize the input multi-variable
    auto volumeOperator = std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
        nonlinearMechanicsOperator->getVolumeOperator() );
    auto inputMultiVariable = std::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(
        volumeOperator->getInputVariable() );
    std::vector<AMP::LinearAlgebra::Variable::shared_ptr> inputVariables;
    std::vector<AMP::Discretization::DOFManager::shared_ptr> inputDOFs;
    for ( size_t i = 0; i < inputMultiVariable->numVariables(); i++ ) {
        inputVariables.push_back( inputMultiVariable->getVariable( i ) );
        if ( i == AMP::Operator::Mechanics::DISPLACEMENT )
            inputDOFs.push_back( displDofMap );
        else if ( i == AMP::Operator::Mechanics::TEMPERATURE )
            inputDOFs.push_back( nodalDofMap );
        else if ( i == AMP::Operator::Mechanics::BURNUP )
            inputDOFs.push_back( nodalDofMap );
        else if ( i == AMP::Operator::Mechanics::OXYGEN_CONCENTRATION )
            inputDOFs.push_back( nodalDofMap );
        else if ( i == AMP::Operator::Mechanics::LHGR )
            inputDOFs.push_back( nodalDofMap );
        else if ( i == AMP::Operator::Mechanics::TOTAL_NUMBER_OF_VARIABLES )
            inputDOFs.push_back( nodalDofMap );
        else
            AMP_ERROR( "Unknown variable" );
    }

    // initialize the output variable
    auto outputVariable = nonlinearThermoMechanicsOperator->getOutputVariable();

    // create solution, rhs, and residual vectors
    auto solVec =
        AMP::LinearAlgebra::MultiVector::create( inputMultiVariable, meshAdapter->getComm() );
    auto rhsVec = AMP::LinearAlgebra::MultiVector::create( outputVariable, meshAdapter->getComm() );
    auto resVec = AMP::LinearAlgebra::MultiVector::create( outputVariable, meshAdapter->getComm() );
    for ( size_t i = 0; i < inputVariables.size(); i++ ) {
        if ( inputVariables[i].get() != nullptr ) {
            solVec->addVector(
                AMP::LinearAlgebra::createVector( inputDOFs[i], inputVariables[i] ) );
            rhsVec->addVector(
                AMP::LinearAlgebra::createVector( inputDOFs[i], inputVariables[i] ) );
            resVec->addVector(
                AMP::LinearAlgebra::createVector( inputDOFs[i], inputVariables[i] ) );
        }
    }

    // set up the shift and scale parameters
    double shift[2];
    double scale[2];
    shift[0] = 0.;
    shift[1] = 0.;
    scale[0] = 1.;
    scale[1] = 1.;
    auto transportModel =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>( thermalTransportModel );
    auto matTh         = transportModel->getMaterial();
    auto thermOperator = std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
        nonlinearThermalOperator->getVolumeOperator() );
    if ( thermOperator->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE ) {
        std::string property = "ThermalConductivity";
        if ( ( matTh->property( property ) )->is_argument( "temperature" ) ) {
            auto range =
                ( matTh->property( property ) )->get_arg_range( "temperature" ); // Compile error
            scale[1] = range[1] - range[0];
            shift[1] = range[0] + 0.001 * scale[1];
            scale[1] *= 0.999;
        }
    }

    //----------------------------------------------------------------------------//
    auto thermVar                = nonlinearThermalOperator->getOutputVariable();
    auto referenceTemperatureVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermVar, true );
    referenceTemperatureVec->setToScalar( 300.0 );
    nonlinearMechanicsVolumeOperator->setReferenceTemperature( referenceTemperatureVec );

    // now construct the linear BVP operator for mechanics
    AMP_INSIST( input_db->keyExists( "testLinearMechanicsOperator" ), "key missing!" );
    auto linearMechanicsOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testLinearMechanicsOperator", input_db, mechanicsMaterialModel ) );

    // now construct the linear BVP operator for thermal
    AMP_INSIST( input_db->keyExists( "testLinearThermalOperator" ), "key missing!" );
    auto linearThermalOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testLinearThermalOperator", input_db, thermalTransportModel ) );

    // create a column operator object for linear thermomechanics
    auto linearThermoMechanicsOperator = std::make_shared<AMP::Operator::ColumnOperator>();
    linearThermoMechanicsOperator->append( linearMechanicsOperator );
    linearThermoMechanicsOperator->append( linearThermalOperator );

    ut->passes( exeName + " : create" );

    // test apply
    std::string msgPrefix = exeName + " : apply";
    auto testOperator     = nonlinearThermoMechanicsOperator;
    applyTests( ut, msgPrefix, testOperator, rhsVec, solVec, resVec, shift, scale, 2 );

    ut->passes( msgPrefix );

    auto resetParams = nonlinearThermoMechanicsOperator->getParameters( "Jacobian", solVec );

    ut->passes( exeName + " : getJacobianParameters" );

    linearThermoMechanicsOperator->reset( resetParams );

    std::cout << "Tested the reset function " << std::endl;

    ut->passes( exeName + " : Linear::reset" );
}

int testNonlinearThermoMechanics( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.emplace_back( "nonlinearBVP-Mechanics-ThermalStrain-Thermal-UO2MSRZC09-1" );

    for ( auto &exeName : exeNames )
        thermoMechanicsTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
