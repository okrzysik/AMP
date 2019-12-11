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
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

#include "../applyTests.h"

#include <iostream>
#include <string>


static void thermalOxygenDiffusionTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );


    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    std::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( database ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr manager     = AMP::Mesh::Mesh::buildMesh( meshParams );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = manager->Subset( "brick" );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear thermal
    AMP_INSIST( input_db->keyExists( "testNonlinearThermalOperator" ), "key missing!" );

    std::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator =
        std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearThermalOperator", input_db ) );
    std::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> nonlinearThermalVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator->getVolumeOperator() );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalMaterialModel =
        nonlinearThermalVolumeOperator->getTransportModel();

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear oxygen diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearOxygenOperator" ), "key missing!" );

    std::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearOxygenOperator =
        std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearOxygenOperator", input_db ) );
    std::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> nonlinearOxygenVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearOxygenOperator->getVolumeOperator() );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> oxygenTransportModel =
        nonlinearOxygenVolumeOperator->getTransportModel();

    std::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> fickOperator =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearOxygenOperator->getVolumeOperator() );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a column operator object for nonlinear thermothermal
    std::shared_ptr<AMP::Operator::OperatorParameters> params;
    std::shared_ptr<AMP::Operator::ColumnOperator> nonlinearThermalOxygenOperator(
        new AMP::Operator::ColumnOperator( params ) );
    nonlinearThermalOxygenOperator->append( nonlinearThermalOperator );
    nonlinearThermalOxygenOperator->append( nonlinearOxygenOperator );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the input multi-variable
    std::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> volumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator->getVolumeOperator() );
    std::shared_ptr<AMP::LinearAlgebra::MultiVariable> inputVariable(
        new AMP::LinearAlgebra::MultiVariable( "inputVariable" ) );
    // inputVariable->add(volumeOperator->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));
    // inputVariable->add(volumeOperator->getInputVariable(AMP::Operator::Diffusion::CONCENTRATION));
    std::shared_ptr<AMP::LinearAlgebra::MultiVariable> tmp =
        std::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(
            volumeOperator->getInputVariable() );
    for ( size_t i = 0; i < tmp->numVariables(); i++ ) {
        if ( tmp->getVariable( i ).get() != nullptr )
            inputVariable->add( tmp->getVariable( i ) );
    }

    // initialize the output multi-variable
    AMP::LinearAlgebra::Variable::shared_ptr outputVariable =
        nonlinearThermalOxygenOperator->getOutputVariable();

    // Create a DOF manager for a nodal vector
    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    AMP::Discretization::DOFManager::shared_ptr nodalDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );

    // create solution, rhs, and residual vectors
    AMP::LinearAlgebra::Vector::shared_ptr solVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, inputVariable );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, outputVariable );
    AMP::LinearAlgebra::Vector::shared_ptr resVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, outputVariable );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // set up the frozen variables for each operator
    // first get defaults
    double defTemp, defConc;
    std::shared_ptr<AMP::Operator::DiffusionTransportModel> thermalTransportModel =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>( thermalMaterialModel );
    std::shared_ptr<AMP::Operator::DiffusionTransportModel> oxyModel =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>( oxygenTransportModel );
    defConc = thermalTransportModel->getDefault( AMP::Operator::Diffusion::CONCENTRATION );
    defTemp = thermalTransportModel->getDefault( AMP::Operator::Diffusion::TEMPERATURE );

    // next get vectors
    AMP::LinearAlgebra::Vector::shared_ptr tempVec =
        solVec->subsetVectorForVariable( inputVariable->getVariable( 0 ) );
    AMP::LinearAlgebra::Vector::shared_ptr concVec =
        solVec->subsetVectorForVariable( inputVariable->getVariable( 1 ) );
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
    AMP::Materials::Material::shared_ptr matt = thermalTransportModel->getMaterial();
    AMP::Materials::Material::shared_ptr mato = oxyModel->getMaterial();
    if ( volumeOperator->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE ) {
        std::string property = "ThermalConductivity";
        if ( ( matt->property( property ) )->is_argument( "temperature" ) ) {
            range = ( matt->property( property ) )->get_arg_range( "temperature" ); // Compile error
            scale[0] = range[1] - range[0];
            shift[0] = range[0] + 0.001 * scale[0];
            scale[0] *= 0.999;
        }
    }
    // the Fick has a principal variable of oxygen
    if ( fickOperator->getPrincipalVariableId() == AMP::Operator::Diffusion::CONCENTRATION ) {
        std::string property = "FickCoefficient";
        if ( ( mato->property( property ) )->is_argument( "concentration" ) ) {
            range =
                ( mato->property( property ) )->get_arg_range( "concentration" ); // Compile error
            scale[1] = range[1] - range[0];
            shift[1] = range[0] + 0.001 * scale[1];
            scale[1] *= 0.999;
        }
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for thermal
    AMP_INSIST( input_db->keyExists( "testLinearThermalOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearThermalOperator", input_db, thermalMaterialModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for oxygen
    AMP_INSIST( input_db->keyExists( "testLinearOxygenOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::LinearBVPOperator> linearOxygenOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearOxygenOperator", input_db, oxygenTransportModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a column operator object for linear thermomechanics
    std::shared_ptr<AMP::Operator::ColumnOperator> linearThermalOxygenOperator(
        new AMP::Operator::ColumnOperator( params ) );
    linearThermalOxygenOperator->append( linearThermalOperator );
    linearThermalOxygenOperator->append( linearOxygenOperator );

    ut->passes( exeName + " : create" );

    // test apply
    std::string msgPrefix = exeName + " : apply";
    std::shared_ptr<AMP::Operator::Operator> testOperator =
        std::dynamic_pointer_cast<AMP::Operator::Operator>( nonlinearThermalOxygenOperator );
    applyTests( ut, msgPrefix, testOperator, rhsVec, solVec, resVec, shift, scale, 2 );

    ut->passes( msgPrefix );

    std::shared_ptr<AMP::Operator::OperatorParameters> resetParams =
        nonlinearThermalOxygenOperator->getParameters( "Jacobian", solVec );

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
