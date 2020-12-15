#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshElement.h"
#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/ampmesh/structured/BoxMesh.h"
#include "AMP/ampmesh/testHelpers/meshTests.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"

#include "meshGenerators.h"

#include "ProfilerApp.h"


template<class GENERATOR>
void runTest( AMP::UnitTest &ut )
{
    auto generator = std::make_shared<GENERATOR>();
    generator->build_mesh();
    AMP::Mesh::meshTests::MeshPerformance( &ut, generator->getMesh() );
}


void testMeshGenerators( AMP::UnitTest &ut )
{
    PROFILE_START( "testMeshGenerators" );
    // AMP mesh generators
    runTest<AMP::unit_test::AMPCubeGenerator<4>>( ut );
    runTest<AMP::unit_test::AMPCylinderGenerator>( ut );
    runTest<AMP::unit_test::AMPMultiMeshGenerator>( ut );
    runTest<AMP::unit_test::AMPCubeGenerator<4>>( ut );
// libMesh generators
#ifdef USE_EXT_LIBMESH
    runTest<AMP::unit_test::LibMeshCubeGenerator<5>>( ut );
    runTest<AMP::unit_test::libMeshThreeElementGenerator>( ut );
#ifdef USE_AMP_DATA
    runTest<AMP::unit_test::ExodusReaderGenerator<>>( ut );
    runTest<AMP::unit_test::ExodusReaderGenerator<2>>( ut );
#endif
#endif
    PROFILE_STOP( "testMeshGenerators" );
}


void testInputMesh( AMP::UnitTest &ut, std::string filename )
{
    PROFILE_START( "testInputMesh" );
    // Read the input file
    auto input_db = AMP::Database::parseInputFile( filename );

    // Get the Mesh database and create the mesh parameters
    auto database = input_db->getDatabase( "Mesh" );
    auto params   = std::make_shared<AMP::Mesh::MeshParameters>( database );
    params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );

    // Create the meshes from the input database
    auto mesh = AMP::Mesh::Mesh::buildMesh( params );

    // Run the mesh tests
    AMP::Mesh::meshTests::MeshPerformance( &ut, mesh );
    ;
    PROFILE_STOP( "testInputMesh" );
}


// Main function
int main( int argc, char **argv )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;
    PROFILE_ENABLE();
    PROFILE_START( "Run tests" );

    if ( argc == 1 ) {
        // Run the default tests
        testMeshGenerators( ut );
    } else {
        // Test each given input file
        for ( int i = 1; i < argc; i++ )
            testInputMesh( ut, argv[i] );
    }

    // Save the timing results
    PROFILE_STOP( "Run tests" );
    PROFILE_SAVE( "test_MeshPerformance" );

    // Print the results and return
    ut.report();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
