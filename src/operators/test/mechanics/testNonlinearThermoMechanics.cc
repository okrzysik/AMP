
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"

#include <iostream>
#include <string>

#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/vectors/VectorBuilder.h"

#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/LinearOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"

#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"

#include "../applyTests.h"


static void thermoMechanicsTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );


    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    std::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    //-----------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear mechanics
    AMP_INSIST( input_db->keyExists( "testNonlinearMechanicsOperator" ), "key missing!" );

    std::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearMechanicsOperator =
        std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearMechanicsOperator", input_db ) );

    std::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> nonlinearMechanicsVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsOperator->getVolumeOperator() );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> mechanicsMaterialModel =
        nonlinearMechanicsVolumeOperator->getMaterialModel();

    //---------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear thermal diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearThermalOperator" ), "key missing!" );

    std::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator =
        std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearThermalOperator", input_db ) );

    std::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> nonlinearThermalVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator->getVolumeOperator() );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel =
        nonlinearThermalVolumeOperator->getTransportModel();

    //--------------------------------------------------------------------------//
    // create a column operator object for nonlinear thermomechanics
    std::shared_ptr<AMP::Operator::OperatorParameters> params;
    std::shared_ptr<AMP::Operator::ColumnOperator> nonlinearThermoMechanicsOperator(
        new AMP::Operator::ColumnOperator( params ) );
    nonlinearThermoMechanicsOperator->append( nonlinearMechanicsOperator );
    nonlinearThermoMechanicsOperator->append( nonlinearThermalOperator );


    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // Create the relavent DOF managers
    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    AMP::Discretization::DOFManager::shared_ptr nodalDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
    int displacementDOFsPerNode = 3;
    AMP::Discretization::DOFManager::shared_ptr displDofMap =
        AMP::Discretization::simpleDOFManager::create( meshAdapter,
                                                       AMP::Mesh::GeomType::Vertex,
                                                       nodalGhostWidth,
                                                       displacementDOFsPerNode,
                                                       split );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the input multi-variable
    std::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> volumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsOperator->getVolumeOperator() );
    std::shared_ptr<AMP::LinearAlgebra::MultiVariable> inputMultiVariable =
        std::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(
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
    AMP::LinearAlgebra::Variable::shared_ptr outputVariable =
        nonlinearThermoMechanicsOperator->getOutputVariable();

    // create solution, rhs, and residual vectors
    std::shared_ptr<AMP::LinearAlgebra::MultiVector> solVec =
        AMP::LinearAlgebra::MultiVector::create( inputMultiVariable, meshAdapter->getComm() );
    std::shared_ptr<AMP::LinearAlgebra::MultiVector> rhsVec =
        AMP::LinearAlgebra::MultiVector::create( outputVariable, meshAdapter->getComm() );
    std::shared_ptr<AMP::LinearAlgebra::MultiVector> resVec =
        AMP::LinearAlgebra::MultiVector::create( outputVariable, meshAdapter->getComm() );
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

    //-----------------------------------------------------------------------//
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
    AMP::LinearAlgebra::Variable::shared_ptr thermVar =
        nonlinearThermalOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr referenceTemperatureVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, thermVar, true );
    referenceTemperatureVec->setToScalar( 300.0 );
    nonlinearMechanicsVolumeOperator->setReferenceTemperature( referenceTemperatureVec );

    //---------------------------------------------------------------------------//
    // now construct the linear BVP operator for mechanics
    AMP_INSIST( input_db->keyExists( "testLinearMechanicsOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::LinearBVPOperator> linearMechanicsOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearMechanicsOperator", input_db, mechanicsMaterialModel ) );

    //--------------------------------------------------------------------------//
    // now construct the linear BVP operator for thermal
    AMP_INSIST( input_db->keyExists( "testLinearThermalOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearThermalOperator", input_db, thermalTransportModel ) );

    //-------------------------------------------------------------------------//
    // create a column operator object for linear thermomechanics
    std::shared_ptr<AMP::Operator::ColumnOperator> linearThermoMechanicsOperator(
        new AMP::Operator::ColumnOperator( params ) );
    linearThermoMechanicsOperator->append( linearMechanicsOperator );
    linearThermoMechanicsOperator->append( linearThermalOperator );

    ut->passes( exeName + " : create" );

    // test apply
    std::string msgPrefix                                 = exeName + " : apply";
    std::shared_ptr<AMP::Operator::Operator> testOperator = nonlinearThermoMechanicsOperator;
    applyTests( ut, msgPrefix, testOperator, rhsVec, solVec, resVec, shift, scale, 2 );

    ut->passes( msgPrefix );

    std::shared_ptr<AMP::Operator::OperatorParameters> resetParams =
        nonlinearThermoMechanicsOperator->getParameters( "Jacobian", solVec );

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
