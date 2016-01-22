#include "utils/AMPManager.h"
#include "materials/Material.h"
#include "operators/NeutronicsRhs.h"
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

#include "operators/subchannel/FlowFrapconOperator.h"
#include "utils/Writer.h"
#include "vectors/Vector.h"

#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "vectors/SimpleVector.h"

#define ITFAILS ut.failure( __LINE__ );
#define UNIT_TEST( a ) \
    if ( !( a ) )      \
        ut.failure( __LINE__ );

void adjust( const AMP::LinearAlgebra::Vector::shared_ptr vec,
             AMP::LinearAlgebra::Vector::shared_ptr work )
{
    work->setToScalar( 301. );
    AMP::LinearAlgebra::Vector &x = *vec;
    AMP::LinearAlgebra::Vector &y = *work;
    vec->add( x, y );
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

void flowTest( AMP::UnitTest *ut )
{
    // Input and output file names
    //  #include <string>
    std::string exeName( "testFlowOperator" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    ////////////////////////////////////
    //    INITIALIZE THE PROBLEM      //
    ////////////////////////////////////


    // Construct a smart pointer to a new database.
    //  #include "utils/shared_ptr.h"
    //  #include "utils/InputDatabase.h"
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );

    // Fill the database from the input file.
    //  #include "utils/InputManager.h"
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    // Print from all cores into the output files
    //   #include "utils/PIO.h"
    AMP::PIO::logAllNodes( log_file );

    // Get the Mesh database and create the mesh parameters
    AMP::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( database ) );
    params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr manager     = AMP::Mesh::Mesh::buildMesh( params );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = manager->Subset( "bar" );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    //////////////////////////////////////
    //     CREATE THE FLOW OPERATOR     //
    //////////////////////////////////////

    AMP_INSIST( input_db->keyExists( "FlowFrapconOperator" ),
                "Key ''FlowFrapconOperator'' is missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> flowtransportModel;
    AMP::shared_ptr<AMP::Operator::FlowFrapconOperator> flowOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::FlowFrapconOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "FlowFrapconOperator", input_db, flowtransportModel ) );

    AMP::LinearAlgebra::Variable::shared_ptr inputVariable  = flowOperator->getInputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr outputVariable = flowOperator->getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr solVec =
        AMP::LinearAlgebra::SimpleVector<double>::create( 10, inputVariable );
    AMP::LinearAlgebra::Vector::shared_ptr cladVec =
        AMP::LinearAlgebra::SimpleVector<double>::create( 10, inputVariable );

    AMP::LinearAlgebra::Vector::shared_ptr rhsVec =
        AMP::LinearAlgebra::SimpleVector<double>::create( 10, outputVariable );
    AMP::LinearAlgebra::Vector::shared_ptr resVec =
        AMP::LinearAlgebra::SimpleVector<double>::create( 10, outputVariable );
    AMP::LinearAlgebra::Vector::shared_ptr workVec =
        AMP::LinearAlgebra::SimpleVector<double>::create( 10, inputVariable );

    cladVec->setToScalar( 300 );
    flowOperator->setVector( cladVec );

    ut->passes( exeName + ": create" );

    AMP::shared_ptr<AMP::Operator::Operator> testOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::Operator>( flowOperator );
    //  flowOperator->apply(nullVec, TemperatureInKelvinVec, flowOutputVec, 1., 0.);

    // test apply
    std::string msgPrefix = exeName + ": apply";
    applyTest( ut, msgPrefix, testOperator, rhsVec, solVec, resVec, workVec );

    ut->passes( msgPrefix );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );

    AMP::UnitTest ut;

    flowTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
