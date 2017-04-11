#include "utils/AMPManager.h"
#include "ampmesh/Mesh.h"
#include "ampmesh/StructuredMeshHelper.h"
#include "ampmesh/structured/BoxMesh.h"
#include "utils/AMP_MPI.h"
#include "utils/UnitTest.h"

#include "meshGenerators.h"


void runTest( AMP::UnitTest *ut )
{
    // Create a simple structured mesh
    size_t N_faces_tot = 227;
    AMP::unit_test::AMPCubeGenerator3<3, 4, 5> generator;
    generator.build_mesh();
    AMP::Mesh::Mesh::shared_ptr mesh = generator.getMesh();

    // Get iterators over each of the face types
    AMP::Mesh::MeshIterator faces   = mesh->getIterator( AMP::Mesh::GeomType::Face, 0 );
    AMP::Mesh::MeshIterator x_faces = AMP::Mesh::StructuredMeshHelper::getYZFaceIterator( mesh, 0 );
    AMP::Mesh::MeshIterator y_faces = AMP::Mesh::StructuredMeshHelper::getXZFaceIterator( mesh, 0 );
    AMP::Mesh::MeshIterator z_faces = AMP::Mesh::StructuredMeshHelper::getXYFaceIterator( mesh, 0 );
    size_t N_faces                  = mesh->getComm().sumReduce( faces.size() );
    size_t N_x_faces                = mesh->getComm().sumReduce( x_faces.size() );
    size_t N_y_faces                = mesh->getComm().sumReduce( y_faces.size() );
    size_t N_z_faces                = mesh->getComm().sumReduce( z_faces.size() );
    if ( N_faces_tot == N_faces )
        ut->passes( "Total number of faces match" );
    else
        ut->failure( "Total number of faces match" );
    if ( N_x_faces + N_y_faces + N_z_faces == N_faces )
        ut->passes( "Number of faces match" );
    else
        ut->failure( "Number of faces match" );

    if ( mesh->getComm().getSize() > 1 ) {
        bool pass = false;
        for ( int d = 0; d < 3; d++ ) {
            AMP::Mesh::MeshIterator face0 =
                AMP::Mesh::StructuredMeshHelper::getYZFaceIterator( mesh, 0 );
            AMP::Mesh::MeshIterator face1 =
                AMP::Mesh::StructuredMeshHelper::getYZFaceIterator( mesh, 1 );
            if ( face1.size() > face0.size() )
                pass = true;
        }
        if ( pass )
            ut->passes( "Found ghost faces on all processors" );
        else
            ut->failure( "Found ghost faces on all processors" );
    }
}


// Main function
int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    runTest( &ut );
    ut.report();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
