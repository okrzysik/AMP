#include "utils/AMPManager.h"
#include "materials/Material.h"
#include "utils/AMPManager.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/Utilities.h"
#include "utils/shared_ptr.h"
#include "vectors/Variable.h"
#include <string>

#include "ProfilerApp.h"
#include "operators/ElementOperationFactory.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/libmesh/SourceNonlinearElement.h"
#include "operators/libmesh/VolumeIntegralOperator.h"
#include "utils/Writer.h"
#include "vectors/Vector.h"

#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "ampmesh/Mesh.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include <exception>


void resetTests( AMP::UnitTest *ut,
                 std::string msgPrefix,
                 AMP::shared_ptr<AMP::Mesh::Mesh>,
                 AMP::shared_ptr<AMP::Operator::Operator>,
                 AMP::shared_ptr<AMP::InputDatabase> )
{
    ut->passes( msgPrefix );
    std::cout.flush();
}

void adjust( const AMP::LinearAlgebra::Vector::shared_ptr vec,
             AMP::LinearAlgebra::Vector::shared_ptr work )
{
    work->setToScalar( 301. );
    AMP::LinearAlgebra::Vector &x = *vec;
    AMP::LinearAlgebra::Vector &y = *work;
    vec->add( x, y );
    vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

void applyTest( AMP::UnitTest *ut,
                std::string msgPrefix,
                AMP::shared_ptr<AMP::Operator::Operator> &testOperator,
                AMP::LinearAlgebra::Vector::shared_ptr rhsVec,
                AMP::LinearAlgebra::Vector::shared_ptr solVec,
                AMP::LinearAlgebra::Vector::shared_ptr resVec,
                AMP::LinearAlgebra::Vector::shared_ptr workVec )
{
    // first test for apply - random values in all three input vectors
    try {
        for ( int j = 0; j < 3; j++ ) {
            solVec->setRandomValues();
            rhsVec->setRandomValues();
            resVec->setRandomValues();
            adjust( solVec, workVec );
            testOperator->residual( rhsVec, solVec, resVec );
        } // end for j
        ut->passes( msgPrefix + " : apply with random f, u, r, a=1, b=-1.0" );
    } catch ( std::exception ) {
        ut->failure( msgPrefix + " : apply with random f, u, r, a=1, b=-1.0" );
    }

    // second test for apply - f NULL, u, r, random values
    try {
        for ( int j = 0; j < 3; j++ ) {
            AMP::LinearAlgebra::Vector::shared_ptr fVec;
            solVec->setRandomValues();
            resVec->setRandomValues();
            adjust( solVec, workVec );
            testOperator->residual( fVec, solVec, resVec );
        } // end for j
        ut->passes( msgPrefix + " : apply with f NULL, random u, r, a=1, b=-1.0" );
    } catch ( std::exception ) {
        ut->failure( msgPrefix + " : apply with f NULL, random u, r, a=1, b=-1.0" );
    }

    // R.S.: u is allowed to be NULL for some operators. For example, operators
    // with an in-place apply. However, this test is not meant to be used with those operators.
    // third test for apply - u NULL, f, r, random values
    try {
        for ( int j = 0; j < 3; j++ ) {
            AMP::LinearAlgebra::Vector::shared_ptr uVec;
            rhsVec->setRandomValues();
            resVec->setRandomValues();
            testOperator->residual( rhsVec, uVec, resVec );
        } // end for j
        ut->failure( msgPrefix +
                     " : apply with u NULL, random values in the vectors f,r, a=1, b=-1.0" );
    } catch ( std::exception ) {
        ut->passes( msgPrefix +
                    " : apply with u NULL, random values in the vectors f,r, a=1, b=-1.0" );
    }

    // fourth test for apply - r NULL, f, u, random values
    try {
        for ( int j = 0; j < 3; j++ ) {
            AMP::LinearAlgebra::Vector::shared_ptr rVec;
            solVec->setRandomValues();
            rhsVec->setRandomValues();
            adjust( solVec, workVec );
            testOperator->residual( rhsVec, solVec, rVec );
        } // end for j
        ut->failure( msgPrefix +
                     " : apply with r NULL, random values in the vectors f,u, a=1, b=-1.0" );
    } catch ( std::exception ) {
        ut->passes( msgPrefix +
                    " : apply with r NULL, random values in the vectors f,u, a=1, b=-1.0" );
    }

    // fifth test for apply - f NULL, u NULL, r, random values
    try {
        for ( int j = 0; j < 3; j++ ) {
            AMP::LinearAlgebra::Vector::shared_ptr fVec;
            AMP::LinearAlgebra::Vector::shared_ptr uVec;
            resVec->setRandomValues();
            testOperator->residual( fVec, uVec, resVec );
        } // end for j
        ut->failure( msgPrefix +
                     " : apply with f NULL, u NULL random values in the vector r, a=1, b=-1.0" );
    } catch ( std::exception ) {
        ut->passes( msgPrefix +
                    " : apply with f NULL, u NULL random values in the vector r, a=1, b=-1.0" );
    }

    // sixth test for apply - u NULL, r NULL, f, random values
    try {
        for ( int j = 0; j < 3; j++ ) {
            AMP::LinearAlgebra::Vector::shared_ptr uVec;
            AMP::LinearAlgebra::Vector::shared_ptr rVec;
            rhsVec->setRandomValues();
            testOperator->residual( rhsVec, uVec, rVec );
        } // end for j
        ut->failure( msgPrefix +
                     " : apply with u NULL, r NULL, random values in the vector f, a=1, b=-1.0" );
    } catch ( std::exception ) {
        ut->passes( msgPrefix +
                    " : apply with u NULL, r NULL, random values in the vector f, a=1, b=-1.0" );
    }

    // seventh test for apply - r NULL, f NULL, u random values
    try {
        for ( int j = 0; j < 3; j++ ) {
            AMP::LinearAlgebra::Vector::shared_ptr rVec;
            AMP::LinearAlgebra::Vector::shared_ptr fVec;
            solVec->setRandomValues();
            adjust( solVec, workVec );
            testOperator->residual( fVec, solVec, rVec );
        } // end for j
        ut->failure( msgPrefix +
                     " : apply with f, r NULL, random values in the vector u, a=1, b=-1.0" );
    } catch ( std::exception ) {
        ut->passes( msgPrefix +
                    " : apply with f, r NULL, random values in the vector u, a=1, b=-1.0" );
    }

    // eighth test for apply - r NULL, f NULL, u NULL
    try {
        for ( int j = 0; j < 3; j++ ) {
            AMP::LinearAlgebra::Vector::shared_ptr rVec;
            AMP::LinearAlgebra::Vector::shared_ptr fVec;
            AMP::LinearAlgebra::Vector::shared_ptr uVec;
            testOperator->residual( fVec, uVec, rVec );
        } // end for j
        ut->failure( msgPrefix + " : apply with f, u, r NULL, a=1, b=-1.0" );
    } catch ( std::exception ) {
        ut->passes( msgPrefix + " : apply with f, u, r NULL, a=1, b=-1.0" );
    }
#if 0
    // ninth test for apply - random values in all three input vectors, a=0, b=1
    rhsVec->setRandomValues();
    resVec->setRandomValues();
    solVec->setRandomValues();
    adjust(solVec, workVec);
    testOperator->apply(rhsVec, solVec, resVec, 0.0, 1.0);
    rhsVec->subtract(rhsVec, resVec);
    if (AMP::Utilities::approx_equal(rhsVec->L2Norm(), 0.0)) {
        ut->passes(msgPrefix + " : apply with random values in the vectors f,u,r, a=0.0, b=1.0");
    } else {
        ut->failure(msgPrefix + " : apply with random values in the vectors f,u,r, a=0.0, b=1.0");
    }

    // tenth test for apply - random values in all three input vectors, a=0, b=-1, to test scaling
    rhsVec->setRandomValues();
    resVec->setRandomValues();
    solVec->setRandomValues();
    adjust(solVec, workVec);
    testOperator->apply(rhsVec, solVec, resVec, 0.0, -1.0);
    rhsVec->add(rhsVec, resVec);
    if (AMP::Utilities::approx_equal(rhsVec->L2Norm(), 0.0)) {
        ut->passes(msgPrefix + " : apply with random values in the vectors f,u,r, a=0.0, b=-1.0 (test scaling of f)");
    } else {
        ut->failure(msgPrefix + " : apply with random values in the vectors f,u,r, a=0.0, b=-1.0 (test scaling of f)");
    }
#endif
}


void sourceTest( AMP::UnitTest *ut, std::string exeName )
{
    // Initialization
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logAllNodes( log_file );

    std::cout << "testing with input file " << input_file << std::endl;
    std::cout.flush();

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );

    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );

    //--------------------------------------------------
    //   CREATE THE VOLUME INTEGRAL OPERATOR -----------
    //--------------------------------------------------

    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
    AMP::shared_ptr<AMP::Database> sourceDatabase =
        input_db->getDatabase( "VolumeIntegralOperator" );
    AMP::pout << "before sourceOp" << std::endl;
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "VolumeIntegralOperator", input_db, transportModel ) );
    AMP::pout << "after sourceOp" << std::endl;
    AMP::LinearAlgebra::Variable::shared_ptr inputVariable  = sourceOperator->getInputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr outputVariable = sourceOperator->getOutputVariable();

    // Create a DOF manager for a gauss point vector
    int DOFsPerElement  = 8;
    int DOFsPerNode     = 1;
    int ghostWidth      = 0;
    int nodalGhostWidth = 1;
    bool split          = true;
    AMP::Discretization::DOFManager::shared_ptr gaussPointDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::Volume, ghostWidth, DOFsPerElement, split );
    AMP::Discretization::DOFManager::shared_ptr nodalDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth, DOFsPerNode, split );

    AMP::LinearAlgebra::Vector::shared_ptr solVec =
        AMP::LinearAlgebra::createVector( gaussPointDofMap, inputVariable, split );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, outputVariable, split );
    AMP::LinearAlgebra::Vector::shared_ptr resVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, outputVariable, split );
    AMP::LinearAlgebra::Vector::shared_ptr workVec =
        AMP::LinearAlgebra::createVector( gaussPointDofMap, inputVariable, split );

    ut->passes( exeName + " : create" );

    std::string msgPrefix = exeName + ": VolumeIntegralOperator reset";
    // test reset
    try {
        AMP::shared_ptr<AMP::Operator::VolumeIntegralOperatorParameters> volumeIntegralParameters(
            new AMP::Operator::VolumeIntegralOperatorParameters( sourceDatabase ) );
        sourceOperator->reset( volumeIntegralParameters );
        ut->passes( msgPrefix + " : pass " );
    } catch ( std::exception ) {
        ut->failure( msgPrefix + " : fail " );
    }

    //  resetTests(ut, msgPrefix, meshAdapter, sourceOperator, sourceDatabase);

    // test apply
    msgPrefix = exeName + ": VolumeIntegralOperator apply";
    AMP::shared_ptr<AMP::Operator::Operator> testOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::Operator>( sourceOperator );
    applyTest( ut, msgPrefix, testOperator, rhsVec, solVec, resVec, workVec );

    ut->passes( msgPrefix );
}


// Input and output file names
int main( int argc, char *argv[] )
{

    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    PROFILE_ENABLE( 5 );
    global_profiler.ignore_timer_errors( true );
    PROFILE_START( "main" );

    AMP::UnitTest ut;

    const int NUMFILES          = 3;
    std::string files[NUMFILES] = { "testVolumeIntegral-1",
                                    "testVolumeIntegral-2",
                                    "testVolumeIntegral-3" };

    for ( auto &file : files )
        sourceTest( &ut, file );

    ut.report();
    PROFILE_STOP( "main" );
    PROFILE_SAVE( "testVolumeIntegral-1" );

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
