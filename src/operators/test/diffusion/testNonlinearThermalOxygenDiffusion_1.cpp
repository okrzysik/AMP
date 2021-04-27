#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshParameters.h"
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
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include "../applyTests.h"

#include <iostream>
#include <memory>
#include <string>


static void thermalOxygenDiffusionTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::logOnlyNodeZero( log_file );

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    auto database   = input_db->getDatabase( "Mesh" );
    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( database );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );

    // Create the meshes from the input database
    auto manager     = AMP::Mesh::Mesh::buildMesh( meshParams );
    auto meshAdapter = manager->Subset( "brick" );

    // create a nonlinear BVP operator for nonlinear thermal
    AMP_INSIST( input_db->keyExists( "testNonlinearThermalOperator" ), "key missing!" );

    auto nonlinearThermalOperator = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testNonlinearThermalOperator", input_db ) );
    auto nonlinearThermalVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator->getVolumeOperator() );
    auto thermalMaterialModel = nonlinearThermalVolumeOperator->getTransportModel();

    // create a nonlinear BVP operator for nonlinear oxygen diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearOxygenOperator" ), "key missing!" );

    auto nonlinearOxygenOperator = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testNonlinearOxygenOperator", input_db ) );
    auto nonlinearOxygenVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearOxygenOperator->getVolumeOperator() );
    auto oxygenTransportModel = nonlinearOxygenVolumeOperator->getTransportModel();

    auto fickOperator = std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
        nonlinearOxygenOperator->getVolumeOperator() );

    // create a column operator object for nonlinear thermothermal
    auto nonlinearThermalOxygenOperator = std::make_shared<AMP::Operator::ColumnOperator>();
    nonlinearThermalOxygenOperator->append( nonlinearThermalOperator );
    nonlinearThermalOxygenOperator->append( nonlinearOxygenOperator );

    // initialize the input multi-variable
    auto volumeOperator = std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
        nonlinearThermalOperator->getVolumeOperator() );
    auto inputVariable = std::make_shared<AMP::LinearAlgebra::MultiVariable>( "inputVariable" );
    // inputVariable->add(volumeOperator->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));
    // inputVariable->add(volumeOperator->getInputVariable(AMP::Operator::Diffusion::CONCENTRATION));
    auto tmp = std::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(
        volumeOperator->getInputVariable() );
    for ( size_t i = 0; i < tmp->numVariables(); i++ ) {
        if ( tmp->getVariable( i ) )
            inputVariable->add( tmp->getVariable( i ) );
    }

    // initialize the output multi-variable
    auto outputVariable = nonlinearThermalOxygenOperator->getOutputVariable();

    // Create a DOF manager for a nodal vector
    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    auto nodalDofMap    = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );

    // create solution, rhs, and residual vectors
    auto solVec = AMP::LinearAlgebra::createVector( nodalDofMap, inputVariable );
    auto rhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, outputVariable );
    auto resVec = AMP::LinearAlgebra::createVector( nodalDofMap, outputVariable );

    // set up the frozen variables for each operator
    // first get defaults
    double defTemp, defConc;
    auto thermalTransportModel =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>( thermalMaterialModel );
    auto oxyModel =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>( oxygenTransportModel );
    defConc = thermalTransportModel->getDefault( AMP::Operator::Diffusion::CONCENTRATION );
    defTemp = thermalTransportModel->getDefault( AMP::Operator::Diffusion::TEMPERATURE );

    // next get vectors
    auto tempVec = solVec->subsetVectorForVariable( inputVariable->getVariable( 0 ) );
    auto concVec = solVec->subsetVectorForVariable( inputVariable->getVariable( 1 ) );
    tempVec->setToScalar( defTemp );
    concVec->setToScalar( defConc );

    // set up the shift and scale parameters
    double shift[2];
    double scale[2];
    shift[0]  = 0.;
    shift[1]  = 0.;
    scale[0]  = 1.;
    scale[1]  = 1.;
    auto matt = thermalTransportModel->getMaterial();
    auto mato = oxyModel->getMaterial();
    if ( volumeOperator->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE ) {
        std::string property = "ThermalConductivity";
        if ( ( matt->property( property ) )->is_argument( "temperature" ) ) {
            auto range =
                ( matt->property( property ) )->get_arg_range( "temperature" ); // Compile error
            scale[0] = range[1] - range[0];
            shift[0] = range[0] + 0.001 * scale[0];
            scale[0] *= 0.999;
        }
    }
    // the Fick has a principal variable of oxygen
    if ( fickOperator->getPrincipalVariableId() == AMP::Operator::Diffusion::CONCENTRATION ) {
        std::string property = "FickCoefficient";
        if ( ( mato->property( property ) )->is_argument( "concentration" ) ) {
            auto range =
                ( mato->property( property ) )->get_arg_range( "concentration" ); // Compile error
            scale[1] = range[1] - range[0];
            shift[1] = range[0] + 0.001 * scale[1];
            scale[1] *= 0.999;
        }
    }

    // now construct the linear BVP operator for thermal
    AMP_INSIST( input_db->keyExists( "testLinearThermalOperator" ), "key missing!" );
    auto linearThermalOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testLinearThermalOperator", input_db, thermalMaterialModel ) );

    // now construct the linear BVP operator for oxygen
    AMP_INSIST( input_db->keyExists( "testLinearOxygenOperator" ), "key missing!" );
    auto linearOxygenOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testLinearOxygenOperator", input_db, oxygenTransportModel ) );

    // create a column operator object for linear thermomechanics
    auto linearThermalOxygenOperator = std::make_shared<AMP::Operator::ColumnOperator>();
    linearThermalOxygenOperator->append( linearThermalOperator );
    linearThermalOxygenOperator->append( linearOxygenOperator );

    ut->passes( exeName + " : create" );

    // test apply
    std::string msgPrefix = exeName + " : apply";
    auto testOperator =
        std::dynamic_pointer_cast<AMP::Operator::Operator>( nonlinearThermalOxygenOperator );
    applyTests( ut, msgPrefix, testOperator, rhsVec, solVec, resVec, shift, scale, 2 );

    ut->passes( msgPrefix );

    auto resetParams = nonlinearThermalOxygenOperator->getParameters( "Jacobian", solVec );

    ut->passes( exeName + " : getJacobianParameters" );

    linearThermalOxygenOperator->reset( resetParams );

    ut->passes( exeName + " : Linear::reset" );
}

int testNonlinearThermalOxygenDiffusion_1( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );

    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.emplace_back( "nonlinearBVP-Thermal-Oxygen-UO2MSRZC09-1" );

    for ( auto &exeName : exeNames )
        thermalOxygenDiffusionTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
