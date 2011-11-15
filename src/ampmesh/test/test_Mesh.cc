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

// Function to test the creation/destruction of empty meshes
void testlibMesh( AMP::UnitTest *ut )
{
    // Create a generic MeshParameters object
    boost::shared_ptr<AMP::MemoryDatabase> database(new AMP::MemoryDatabase("Mesh"));
    database->putInteger("dim",3);
    database->putString("FileName","pellet_1x.e");
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
    params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));

    // Create an libMesh mesh
    boost::shared_ptr<AMP::Mesh::libMesh> mesh1(new AMP::Mesh::libMesh(params));    

    // Run the mesh tests
    MeshTestLoop( ut, mesh1 );

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

}



// Main function
int main ( int argc , char ** argv )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = true;
    AMP::AMPManager::startup(argc,argv,startup_properties);
    AMP::UnitTest ut;

    // Run tests on a libmesh mesh
    //testlibMesh( &ut );

    // Run tests on the input file
    testInputMesh( &ut, "input_Mesh" );

    ut.report ();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
