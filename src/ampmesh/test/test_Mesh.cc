#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/AMP_MPI.h"
#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "ampmesh/libmesh/libMesh.h"
#include "utils/MemoryDatabase.h"
#include "meshTestLoop.h"

// Function to test the creation/destruction of empty meshes
void testMeshCreation( AMP::UnitTest *ut )
{
    // Create a generic MeshParameters object
    boost::shared_ptr<AMP::MemoryDatabase> database(new AMP::MemoryDatabase("Mesh"));
    database->putInteger("dim",3);
    database->putString("Filename","pellet_1x.e");
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
    params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));

    // Create an libMesh mesh
    boost::shared_ptr<AMP::Mesh::libMesh> mesh1(new AMP::Mesh::libMesh(params));    

    // Run the mesh tests
    MeshTestLoop( ut, mesh1 );

}



// Main function
int main ( int argc , char ** argv )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);
    AMP::UnitTest ut;

    testMeshCreation( &ut );

    ut.report ();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
