#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include <iostream>
#include <string>

#include "utils/shared_ptr.h"

#include "materials/Material.h"
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/MultiVector.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "utils/Writer.h"

#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"

#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"

#include "operators/BVPOperatorParameters.h"
#include "operators/ColumnOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "operators/boundary/DirichletVectorCorrection.h"


void myTest( AMP::UnitTest *ut, std::string exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

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

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> mechanicsMaterialModel;
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearMechanicsOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearMechanicsOperator", input_db, mechanicsMaterialModel ) );


    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear thermal diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearThermalOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearThermalOperator", input_db, thermalTransportModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
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
    AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> mechanicsVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsOperator->getVolumeOperator() );
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator->getVolumeOperator() );

    // initialize the input multi-variable
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> inputMultiVariable =
        AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(
            mechanicsVolumeOperator->getInputVariable() );
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
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> solVec =
        AMP::LinearAlgebra::MultiVector::create( inputMultiVariable, meshAdapter->getComm() );
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> resVec =
        AMP::LinearAlgebra::MultiVector::create( outputVariable, meshAdapter->getComm() );
    for ( size_t i = 0; i < inputVariables.size(); i++ ) {
        if ( inputVariables[i].get() != nullptr ) {
            AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>( solVec )->addVector(
                AMP::LinearAlgebra::createVector( inputDOFs[i], inputVariables[i] ) );
            AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>( resVec )->addVector(
                AMP::LinearAlgebra::createVector( inputDOFs[i], inputVariables[i] ) );
        }
    }

    // create the following shared pointers for ease of use
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // IMPORTANT:: call init before proceeding any further on the nonlinear mechanics operator
    //  previous line:
    //  AMP::LinearAlgebra::Vector::shared_ptr referenceTemperatureVec = meshAdapter->createVector(
    //  thermalVolumeOperator->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE) );
    //  converted line:
    //  AMP::LinearAlgebra::Vector::shared_ptr referenceTemperatureVec =
    //  AMP::LinearAlgebra::createVector( nodalDofMap,
    //  (thermalVolumeOperator->getInputVariable())->getVariable(AMP::Operator::Diffusion::TEMPERATURE)
    //  );
    //  placeholder line till diffusion is converted:
    AMP::LinearAlgebra::Variable::shared_ptr temperatureVariable =
        AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(
            thermalVolumeOperator->getInputVariable() )
            ->getVariable( AMP::Operator::Diffusion::TEMPERATURE );
    AMP::LinearAlgebra::Vector::shared_ptr referenceTemperatureVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, temperatureVariable );
    referenceTemperatureVec->setToScalar( 300.0 );
    mechanicsVolumeOperator->setReferenceTemperature( referenceTemperatureVec );

    nonlinearMechanicsOperator->apply( solVec, resVec );

    ut->passes( exeName );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back( "testNonlinearThermoMechanics-bug-1" );

    for ( auto &exeName : exeNames ) {
        try {
            myTest( &ut, exeName );
        } catch ( std::exception &err ) {
            std::cout << "ERROR: While testing " << argv[0] << err.what() << std::endl;
            ut.failure( "ERROR: While testing" );
        } catch ( ... ) {
            std::cout << "ERROR: While testing " << argv[0] << "An unknown exception was thrown."
                      << std::endl;
            ut.failure( "ERROR: While testing" );
        }
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
