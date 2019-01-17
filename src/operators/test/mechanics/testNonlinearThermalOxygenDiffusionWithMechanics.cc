#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
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
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include "../applyTests.h"

#include <iostream>
#include <string>


static void thermoMechanicsTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logAllNodes( log_file );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    const int rank = AMP::AMP_MPI( AMP_COMM_WORLD ).getRank();

    //--------------------------------------------------
    //   Create the Mesh.
    //--------------------------------------------------
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
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

    //----------------------------------------------------------------------------------------------------------------------------------------------//
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

    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator->getVolumeOperator() );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear oxygen diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearOxygenOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearOxygenOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearOxygenOperator", input_db ) );
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> nonlinearOxygenVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearOxygenOperator->getVolumeOperator() );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> oxygenTransportModel =
        nonlinearOxygenVolumeOperator->getTransportModel();

    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> fickOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearOxygenOperator->getVolumeOperator() );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a column operator object for nonlinear thermal and oxygen diffusion, and mechanics
    AMP::shared_ptr<AMP::Operator::OperatorParameters> params;
    AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearThermalOxygenDiffusionMechanicsOperator(
        new AMP::Operator::ColumnOperator( params ) );
    nonlinearThermalOxygenDiffusionMechanicsOperator->append( nonlinearMechanicsOperator );
    nonlinearThermalOxygenDiffusionMechanicsOperator->append( nonlinearThermalOperator );
    nonlinearThermalOxygenDiffusionMechanicsOperator->append( nonlinearOxygenOperator );

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
        nonlinearThermalOxygenDiffusionMechanicsOperator->getOutputVariable();

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


    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // set up the frozen variables for each operator
    // first get defaults
    double defTemp, defConc;
    AMP::shared_ptr<AMP::Operator::DiffusionTransportModel> transportModelTh =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>( thermalTransportModel );
    defConc = transportModelTh->getDefault( AMP::Operator::Diffusion::CONCENTRATION );
    AMP::shared_ptr<AMP::Operator::DiffusionTransportModel> transportModelOx =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>( oxygenTransportModel );
    defTemp = transportModelOx->getDefault( AMP::Operator::Diffusion::TEMPERATURE );
    // next get vectors
    AMP::LinearAlgebra::Vector::shared_ptr tempVec =
        solVec->subsetVectorForVariable( inputVariables[AMP::Operator::Mechanics::TEMPERATURE] );
    AMP::LinearAlgebra::Vector::shared_ptr concVec = solVec->subsetVectorForVariable(
        inputVariables[AMP::Operator::Mechanics::OXYGEN_CONCENTRATION] );
    tempVec->setToScalar( defTemp );
    concVec->setToScalar( defConc );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // set up the shift and scale parameters
    double shift[2];
    double scale[2];
    shift[0] = 0.;
    shift[1] = 0.;
    scale[0] = 1.;
    scale[1] = 1.;
    std::vector<double> range( 2 );
    AMP::Materials::Material::shared_ptr matTh = transportModelTh->getMaterial();
    AMP::Materials::Material::shared_ptr matOx = transportModelOx->getMaterial();
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
    // the Fick has a principal variable of oxygen
    if ( fickOperator->getPrincipalVariableId() == AMP::Operator::Diffusion::CONCENTRATION ) {
        std::string property = "FickCoefficient";
        if ( ( matOx->property( property ) )->is_argument( "concentration" ) ) {
            range =
                ( matOx->property( property ) )->get_arg_range( "concentration" ); // Compile error
            scale[0] = range[1] - range[0];
            shift[0] = range[0] + 0.001 * scale[0];
            scale[0] *= 0.999;
        }
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // IMPORTANT:: call init before proceeding any further on the nonlinear mechanics operator
    AMP::LinearAlgebra::Vector::shared_ptr referenceTemperatureVec =
        AMP::LinearAlgebra::createVector(
            nodalDofMap, inputMultiVariable->getVariable( AMP::Operator::Mechanics::TEMPERATURE ) );
    referenceTemperatureVec->setToScalar( 300.0 );
    volumeOperator->setReferenceTemperature( referenceTemperatureVec );
    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for mechanics
    AMP_INSIST( input_db->keyExists( "testLinearMechanicsOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearMechanicsOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearMechanicsOperator", input_db, mechanicsMaterialModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for thermal
    AMP_INSIST( input_db->keyExists( "testLinearThermalOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearThermalOperator", input_db, thermalTransportModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for oxygen
    AMP_INSIST( input_db->keyExists( "testLinearOxygenOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearOxygenOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearOxygenOperator", input_db, oxygenTransportModel ) );
    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a column operator object for linear thermomechanics
    AMP::shared_ptr<AMP::Operator::ColumnOperator> linearThermalOxygenDiffusionMechanicsOperator(
        new AMP::Operator::ColumnOperator( params ) );
    linearThermalOxygenDiffusionMechanicsOperator->append( linearMechanicsOperator );
    linearThermalOxygenDiffusionMechanicsOperator->append( linearThermalOperator );
    linearThermalOxygenDiffusionMechanicsOperator->append( linearOxygenOperator );

    ut->passes( exeName + " : create" );

    // test apply
    std::string msgPrefix = exeName + " : apply";
    AMP::shared_ptr<AMP::Operator::Operator> testOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::Operator>(
            nonlinearThermalOxygenDiffusionMechanicsOperator );
    if ( rank == 0 )
        std::cout << "Running apply tests" << std::endl;
    applyTests( ut, msgPrefix, testOperator, rhsVec, solVec, resVec, shift, scale, 3 );
    AMP::AMP_MPI( AMP_COMM_WORLD ).barrier();
    if ( rank == 0 )
        std::cout << "Finished apply tests" << std::endl;

    ut->passes( msgPrefix );

    AMP::shared_ptr<AMP::Operator::OperatorParameters> resetParams =
        nonlinearThermalOxygenDiffusionMechanicsOperator->getParameters( "Jacobian", solVec );

    ut->passes( exeName + " : getJacobianParameters" );

    linearThermalOxygenDiffusionMechanicsOperator->reset( resetParams );

    ut->passes( exeName + " : Linear::reset" );

    AMP::AMP_MPI( AMP_COMM_WORLD ).barrier();
    if ( rank == 0 )
        std::cout << "Finished tests: " << exeName << std::endl;
}


int testNonlinearThermalOxygenDiffusionWithMechanics( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;
    ut.verbose();

    std::vector<std::string> exeNames;
    exeNames.emplace_back( "nonlinearBVP-Mechanics-ThermalStrain-Thermal-Oxygen-UO2MSRZC09-1" );
    // exeNames.push_back("testNonlinearMechanics-1-reduced");

    for ( auto &exeName : exeNames ) {
        thermoMechanicsTest( &ut, exeName );
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
