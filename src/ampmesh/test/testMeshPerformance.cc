#include "meshGenerators.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "ampmesh/structured/BoxMesh.h"
#include "ampmesh/testHelpers/meshTests.h"

#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/MemoryDatabase.h"
#include "utils/UnitTest.h"


#include "ProfilerApp.h"


template<class GENERATOR>
void runTest( AMP::UnitTest &ut )
{
    auto generator = AMP::make_shared<GENERATOR>();
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
// libmesh generators
#ifdef USE_EXT_LIBMESH
    // Test the libmesh cube generator
    runTest<AMP::unit_test::LibMeshCubeGenerator<5>>( ut );
    runTest<AMP::unit_test::ExodusReaderGenerator<>>( ut );
    runTest<AMP::unit_test::ExodusReaderGenerator<2>>( ut );
    runTest<AMP::unit_test::libMeshThreeElementGenerator>( ut );
#endif
    PROFILE_STOP( "testMeshGenerators" );
}


void testInputMesh( AMP::UnitTest &ut, std::string filename )
{
    PROFILE_START( "testInputMesh" );
    // Read the input file
    auto input_db = AMP::make_shared<AMP::InputDatabase>( "input_db" );
    AMP::InputManager::getManager()->parseInputFile( filename, input_db );

    // Get the Mesh database and create the mesh parameters
    auto database = input_db->getDatabase( "Mesh" );
    auto params   = AMP::make_shared<AMP::Mesh::MeshParameters>( database );
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
