#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/ElementOperationFactory.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/libmesh/RobinMatrixCorrection.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <memory>
#include <string>


static void bcTests( AMP::UnitTest *ut,
                     std::string msgPrefix,
                     std::shared_ptr<AMP::Operator::Operator> &feOperator,
                     std::shared_ptr<AMP::Operator::Operator> &bcOperator,
                     std::shared_ptr<AMP::Database> bcDatabase,
                     AMP::LinearAlgebra::Vector::shared_ptr bcCorrectionVec )
//             std::shared_ptr<AMP::Operator::OperatorParameters> &bcParameters)
{

    try {
        auto tmp_db = std::make_shared<AMP::Database>( "Dummy" );
        tmp_db->putScalar( "skip_params", false );
        tmp_db->putScalar( "print_info_level", 3 );
        tmp_db->putScalar( "alpha", 0.0 );

        auto dummyParameters =
            std::make_shared<AMP::Operator::RobinMatrixCorrectionParameters>( tmp_db );

        bcOperator->reset( dummyParameters );

        ut->failure( "Test" );
    } catch ( const std::exception & ) {

        ut->passes( msgPrefix + ": catches when prefactor alpha is set to zero " );
    }

    ut->passes( msgPrefix );

    try {
        auto bcParameters =
            std::make_shared<AMP::Operator::RobinMatrixCorrectionParameters>( bcDatabase );
        bcParameters->d_inputMatrix =
            ( std::dynamic_pointer_cast<AMP::Operator::LinearFEOperator>( feOperator ) )
                ->getMatrix();
        bcParameters->d_variable = feOperator->getOutputVariable();
        bcOperator->reset( bcParameters );

        bcCorrectionVec->setToScalar( 0.0 );
        ( std::dynamic_pointer_cast<AMP::Operator::BoundaryOperator>( bcOperator ) )
            ->addRHScorrection( bcCorrectionVec );
        AMP_INSIST( ( ( bcCorrectionVec.get() ) != nullptr ), "NULL rhs correction vector" );

        ut->passes( msgPrefix + ": Robin returns a rhs correction vector " );

        // ut.failure(msgPrefix+": BoundaryOperators have changed and this needs to be updated.");
    } catch ( ... ) {

        ut->failure( "Exception" );
    }

    ut->passes( msgPrefix );
    std::cout.flush();
}

static void linearRobinTest( AMP::UnitTest *ut, const std::string &exeName )
{
    // Initialization
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::logAllNodes( log_file );

    std::cout << "testing with input file " << input_file << std::endl;
    std::cout.flush();


    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Get the mesh name
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    std::string mesh_file = input_db->getString( "Mesh" );

    // Create the mesh parameter object
    std::shared_ptr<AMP::Database> database( new AMP::Database( "Mesh" ) );
    database->putScalar( "dim", 3 );
    database->putScalar( "MeshName", "mesh" );
    database->putScalar( "MeshType", "libMesh" );
    database->putScalar( "FileName", mesh_file );
    std::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( database ) );
    params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );

    // Create the mesh
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( params );

    /////////////////////////////////////////////////
    //   CREATE THE LINEAR DIFFUSION BVP OPERATOR  //
    /////////////////////////////////////////////////
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> diffusionOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "DiffusionBVPOperator", input_db, elementModel ) );

    std::shared_ptr<AMP::Database> bcDatabase = std::dynamic_pointer_cast<AMP::Database>(
        input_db->getDatabase( "RobinMatrixCorrection" ) );

    AMP::Operator::Operator::shared_ptr boundaryOp, volumeOp;
    boundaryOp = diffusionOperator->getBoundaryOperator();
    volumeOp   = diffusionOperator->getVolumeOperator();

    // AMP::LinearAlgebra::Vector::shared_ptr bndVec = meshAdapter->createVector(
    // volumeOp->getOutputVariable() );
    std::shared_ptr<AMP::Discretization::DOFManager> NodalScalarDOF =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    AMP::LinearAlgebra::Vector::shared_ptr bndVec =
        AMP::LinearAlgebra::createVector( NodalScalarDOF, volumeOp->getOutputVariable(), true );

    std::string msgPrefix = exeName + "- Boundary Conditions";
    // Test Robin Boundary Conditions
    {
        bcTests( ut, msgPrefix, volumeOp, boundaryOp, bcDatabase, bndVec );
    }

    std::cout.flush();
}
// Input and output file names
int testLinearRobin( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );

    AMP::UnitTest ut;

    const int NUMFILES          = 1;
    std::string files[NUMFILES] = { "LinearOp-Robin-1" };

    for ( auto &file : files )
        linearRobinTest( &ut, file );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
