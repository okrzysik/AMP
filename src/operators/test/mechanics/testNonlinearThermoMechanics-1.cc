
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include <iostream>
#include <string>

#include "discretization/MultiDOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "operators/ColumnOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/LinearOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"

#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"

#include "../applyTests.h"


void thermoMechanicsTest( AMP::UnitTest *ut, std::string exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    //-----------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear mechanics
    AMP_INSIST( input_db->keyExists( "testNonlinearMechanicsOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearMechanicsOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearMechanicsOperator", input_db ) );

    AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> nonlinearMechanicsVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsOperator->getVolumeOperator() );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> mechanicsMaterialModel =
        nonlinearMechanicsVolumeOperator->getMaterialModel();

    //---------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear thermal diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearThermalOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearThermalOperator", input_db ) );

    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> nonlinearThermalVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator->getVolumeOperator() );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel =
        nonlinearThermalVolumeOperator->getTransportModel();

    //--------------------------------------------------------------------------//
    // create a column operator object for nonlinear thermomechanics
    AMP::shared_ptr<AMP::Operator::OperatorParameters> params;
    AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearThermoMechanicsOperator(
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
            meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth, DOFsPerNode, split );
    int displacementDOFsPerNode = 3;
    AMP::Discretization::DOFManager::shared_ptr displDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth, displacementDOFsPerNode, split );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the input multi-variable
    AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> volumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsOperator->getVolumeOperator() );
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> inputMultiVariable =
        AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(
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
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> solVec =
        AMP::LinearAlgebra::MultiVector::create( inputMultiVariable, meshAdapter->getComm() );
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> rhsVec =
        AMP::LinearAlgebra::MultiVector::create( outputVariable, meshAdapter->getComm() );
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> resVec =
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
    std::vector<double> range( 2 );
    AMP::shared_ptr<AMP::Operator::DiffusionTransportModel> transportModel =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>( thermalTransportModel );
    AMP::Materials::Material::shared_ptr matTh = transportModel->getMaterial();
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator->getVolumeOperator() );
    if ( thermOperator->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE ) {
        std::string property = "ThermalConductivity";
        if ( ( matTh->property( property ) )->is_argument( "temperature" ) ) {
            range =
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
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearMechanicsOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearMechanicsOperator", input_db, mechanicsMaterialModel ) );

    //--------------------------------------------------------------------------//
    // now construct the linear BVP operator for thermal
    AMP_INSIST( input_db->keyExists( "testLinearThermalOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearThermalOperator", input_db, thermalTransportModel ) );

    //-------------------------------------------------------------------------//
    // create a column operator object for linear thermomechanics
    AMP::shared_ptr<AMP::Operator::ColumnOperator> linearThermoMechanicsOperator(
        new AMP::Operator::ColumnOperator( params ) );
    linearThermoMechanicsOperator->append( linearMechanicsOperator );
    linearThermoMechanicsOperator->append( linearThermalOperator );

    ut->passes( exeName + " : create" );

    // test apply
    std::string msgPrefix                                 = exeName + " : apply";
    AMP::shared_ptr<AMP::Operator::Operator> testOperator = nonlinearThermoMechanicsOperator;
    applyTests( ut, msgPrefix, testOperator, rhsVec, solVec, resVec, shift, scale, 2 );

    ut->passes( msgPrefix );

    AMP::shared_ptr<AMP::Operator::OperatorParameters> resetParams =
        nonlinearThermoMechanicsOperator->getParameters( "Jacobian", solVec );

    ut->passes( exeName + " : getJacobianParameters" );

    linearThermoMechanicsOperator->reset( resetParams );

    std::cout << "Tested the reset function " << std::endl;

    ut->passes( exeName + " : Linear::reset" );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back( "nonlinearBVP-Mechanics-ThermalStrain-Thermal-UO2MSRZC09-1" );

    for ( auto &exeName : exeNames )
        thermoMechanicsTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
