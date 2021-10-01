#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/ampmesh/structured/PureLogicalMesh.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"

#include "ProfilerApp.h"

#include <cmath>
#include <iomanip>
#include <iostream>


// Main function
int run( const std::array<int, 3> &size, int N_proc )
{

    // Create the mesh database
    auto db = std::make_shared<AMP::Database>( "mesh" );
    db->putScalar( "MeshName", "logical" );
    db->putScalar( "MeshType", "AMP" );
    db->putScalar( "Generator", "logical" );
    db->putVector( "Size", std::vector<int>( size.begin(), size.end() ) );
    db->putScalar( "commSize", N_proc );
    db->putScalar( "commRank", 0 );
    auto params = std::make_shared<AMP::Mesh::MeshParameters>( db );
    params->setComm( AMP_COMM_SELF );
    // Loop through the ranks
    size_t N_node = 0;
    size_t N_elem = 0;
    AMP::Mesh::BoxMesh::Box box;
    for ( int rank = 0; rank < N_proc; rank++ ) {
        db->putScalar( "commRank", rank, AMP::Units(), AMP::Database::Check::Overwrite );
        // Create the mesh
        auto mesh0 = AMP::Mesh::Mesh::buildMesh( params );
        auto mesh  = std::dynamic_pointer_cast<AMP::Mesh::PureLogicalMesh>( mesh0 );
        AMP_ASSERT( mesh );
        // Get the box and check the global size
        box = mesh->getGlobalBox();
        AMP_ASSERT( mesh->numGlobalElements( AMP::Mesh::GeomType::Volume ) == box.size().length() );
        AMP_ASSERT( mesh->numGlobalElements( AMP::Mesh::GeomType::Vertex ) ==
                    ( box.size() + 1 ).length() );
        // Check if we have any empty ranks (will be an error in the future)
        if ( mesh->numLocalElements( AMP::Mesh::GeomType::Volume ) == 0 )
            printf( "Proc %i is empty for %i processors: (%i,%i,%i)\n",
                    rank + 1,
                    N_proc,
                    size[0],
                    size[1],
                    size[2] );
        // Loop through the nodes to make sure we can correctly convert all node indices
        for ( int k = 0; k <= box.last[2] + 1; k++ ) {
            for ( int j = 0; j <= box.last[1] + 1; j++ ) {
                for ( int i = 0; i <= box.last[0] + 1; i++ ) {
                    AMP::Mesh::BoxMesh::MeshElementIndex id1(
                        AMP::Mesh::GeomType::Vertex, 0, i, j, k );
                    auto id2 = mesh->convert( id1 );
                    auto id3 = mesh->convert( id2 );
                    AMP_ASSERT( id1 == id3 );
                    if ( id2.is_local() )
                        N_node++;
                }
            }
        }
        // Loop through the elements to make sure we can correctly convert all node indices
        for ( int k = 0; k <= box.last[2]; k++ ) {
            for ( int j = 0; j <= box.last[1]; j++ ) {
                for ( int i = 0; i <= box.last[0]; i++ ) {
                    AMP::Mesh::BoxMesh::MeshElementIndex id1(
                        AMP::Mesh::GeomType::Volume, 0, i, j, k );
                    auto id2 = mesh->convert( id1 );
                    auto id3 = mesh->convert( id2 );
                    AMP_ASSERT( id1 == id3 );
                    if ( id2.is_local() )
                        N_elem++;
                }
            }
        }
    }
    AMP_ASSERT( N_elem == box.size().length() );
    AMP_ASSERT( N_node == ( box.size() + 1 ).length() );
    return 0;
}


// Main function
int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );

    // Define some problem sizes and number of processors
    std::array<int, 3> size[] = { { 88, 1, 1 }, { 128, 1, 1 }, { 2500, 1, 1 }, { 7, 13, 17 } };
    int N_procs[]             = { 2, 3, 4, 5, 8, 12, 13, 16, 64, 87 };

    // Run the problem
    int N_errors = 0;
    for ( auto s : size ) {
        for ( auto N : N_procs ) {
            N_errors += run( s, N );
        }
    }
    if ( N_errors == 0 )
        std::cout << "All tests passed\n";

    // Shutdown AMP
    AMP::AMPManager::shutdown();
    return N_errors;
}
