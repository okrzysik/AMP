#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/AMP_MPI.h"
#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "ampmesh/libmesh/libMesh.h"
#include "utils/MemoryDatabase.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "meshTestLoop.h"
#include "meshTests.h"
#include "meshGenerators.h"


// Function to test the creation/destruction of a mesh with the mesh generators
// Note: this only runs the mesh tests, not the vector or matrix tests
void testMeshGenerators( AMP::UnitTest *ut )
{
    boost::shared_ptr<AMP::unit_test::MeshGenerator> generator;
    // Test the libmesh cube generator
    generator = boost::shared_ptr<AMP::unit_test::MeshGenerator>( new AMP::unit_test::LibMeshCubeGenerator<5> );
    generator->build_mesh();
    MeshTestLoop( ut, generator->getMesh() );
    // Test the libmesh reader generator
    generator = boost::shared_ptr<AMP::unit_test::MeshGenerator>( new AMP::unit_test::ExodusReaderGenerator);
    generator->build_mesh();
    MeshTestLoop( ut, generator->getMesh() );
    // Test the multimesh generator
    generator = boost::shared_ptr<AMP::unit_test::MeshGenerator>( new AMP::unit_test::MultiMeshGenerator );
    generator->build_mesh();
    MeshTestLoop( ut, generator->getMesh() );
    MeshVectorTestLoop( ut, generator->getMesh() );
    // Test the ThreeElementLGenerator generator
    generator = boost::shared_ptr<AMP::unit_test::MeshGenerator>( new AMP::unit_test::libMeshThreeElementGenerator );
    generator->build_mesh();
    MeshTestLoop( ut, generator->getMesh() );

}


// Function to test the creation/destruction of a libmesh mesh
void testlibMesh( AMP::UnitTest *ut )
{
    // Create a generic MeshParameters object
    boost::shared_ptr<AMP::MemoryDatabase> database(new AMP::MemoryDatabase("Mesh"));
    database->putInteger("dim",3);
    database->putString("MeshName","mesh1");
    database->putString("FileName","pellet_lo_res.e");
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
    params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));

    // Create an libMesh mesh
    boost::shared_ptr<AMP::Mesh::libMesh> mesh(new AMP::Mesh::libMesh(params));    

    // Run the mesh tests
    MeshTestLoop( ut, mesh );
    MeshVectorTestLoop( ut, mesh );
    MeshMatrixTestLoop( ut, mesh );

}


void testInputMesh( AMP::UnitTest *ut, std::string filename )
{
    // Read the input file
    boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( filename , input_db );

    // Get the Mesh database and create the mesh parameters
    boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
    params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));

    // Create the meshes from the input database
    boost::shared_ptr<AMP::Mesh::Mesh> mesh = AMP::Mesh::Mesh::buildMesh(params);

    // Run the mesh tests
    MeshTestLoop( ut, mesh );
    MeshVectorTestLoop( ut, mesh );
    MeshMatrixTestLoop( ut, mesh );

}


void testSubsetMesh( AMP::UnitTest *ut )
{
    // Create the mesh
    boost::shared_ptr<AMP::unit_test::MeshGenerator> generator( 
        new AMP::unit_test::SurfaceSubsetGenerator<AMP::unit_test::ExodusReaderGenerator> );
    generator->build_mesh();
    AMP::Mesh::Mesh::shared_ptr mesh = generator->getMesh();

    // Run the mesh tests
    MeshTestLoop( ut, mesh );
    //MeshVectorTestLoop( ut, mesh );
    //MeshMatrixTestLoop( ut, mesh );

}



// Main function
int main ( int argc , char ** argv )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);
    AMP::UnitTest ut;

    // Run the ID test
    testID( &ut );

    // Run tests on a libmesh mesh
    testlibMesh( &ut );

    // Run tests on the input file
    testInputMesh( &ut, "input_Mesh" );

    // Run the basic tests on all mesh generators
    testMeshGenerators( &ut );

    // Run the tests on the subset meshes
    testSubsetMesh( &ut );

    // Print the results and return
    ut.report ();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
