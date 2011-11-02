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
    /*AMP::Mesh::libMesh test1 = *mesh1;
    AMP::Mesh::Mesh test2 = test1;
    AMP::Mesh::GeomType type = mesh1->getGeomType();
    size_t N1 = test1.numLocalElements(type);
    size_t N2 = test2.numLocalElements(type);
    if ( N1==N2)
        ut->passes("Basic derived class copy works");
    else
        ut->failure("Basic derived class copy works");*/
    

    // Run the mesh tests
    MeshTestLoop( ut, mesh1 );
    
}



// Main function
int main ( int argc , char ** argv )
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    testMeshCreation( &ut );

    ut.report ();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
