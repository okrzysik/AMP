#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/mesh/structured/PureLogicalMesh.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"

#include "ProfilerApp.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>


// Helper function to return the unique elements in an array
template<class TYPE>
static inline std::vector<TYPE> create( const std::string &str, int P = 0 )
{
    AMP::Array<TYPE> x( str );
    std::vector<TYPE> y( x.begin(), x.end() );
    if ( P > 0 ) {
        auto p = AMP::Utilities::primes( P );
        for ( TYPE v : p )
            y.push_back( v );
    }
    AMP::Utilities::unique( y );
    return y;
}


// Function to convert from global index to local index and processor
static inline std::tuple<int, int> convert1( int i, int Np, double n )
{
    int p = std::min<int>( i / n, Np - 1 );
    int k = i - static_cast<int>( p * n );
    return std::tie( p, k );
}


// Function to convert from local index and processor to global index
static inline int convert2( int p, int i, double n ) { return i + static_cast<int>( p * n ); }


// Function to test converting back and forth between indicies
bool testConversion()
{
    bool pass = true;
    auto size = create<int>( "[ 1:500 501:4:1000 10000 100000]", 1000 );
    AMP::Array<int> N_procs( "[ 1:100 127 128 192 256 512 1000 1024 2500 5000 10000 100000 ]" );
    for ( auto x : size ) {
        for ( auto Np : N_procs ) {
            // Compute the average number of elements/processor
            double n = static_cast<double>( x ) / static_cast<double>( Np ) + 1e-12;
            n        = std::max( n, 1.0 );
            for ( int i = 0; i <= x; i++ ) {
                auto [p, i2] = convert1( i, Np, n );
                auto i3      = convert2( p, i2, n );
                if ( i3 != i ) {
                    printf(
                        "Failed conversion (%i,%i,%i) -> (%i,%i) -> %i\n", i, Np, x, p, i2, i3 );
                    pass = false;
                }
            }
        }
    }
    return pass;
}


// Main function
bool run( const std::array<int, 3> &size, const std::array<bool, 3> &periodic, int N_proc )
{
    // Create the mesh database
    auto db = std::make_shared<AMP::Database>( "mesh" );
    db->putScalar( "MeshName", "logical" );
    db->putScalar( "MeshType", "AMP" );
    db->putScalar( "Generator", "logical" );
    db->putVector( "Size", std::vector<int>( size.begin(), size.end() ) );
    db->putVector( "Periodic", std::vector<bool>( periodic.begin(), periodic.end() ) );
    db->putScalar( "commSize", N_proc );
    db->putScalar( "commRank", 0 );
    auto params = std::make_shared<AMP::Mesh::MeshParameters>( db );
    params->setComm( AMP_COMM_SELF );
    // Loop through the ranks
    size_t N_node = 0;
    size_t N_elem = 0;
    AMP::Mesh::BoxMesh::Box globalElemBox, globalNodeBox;
    for ( int rank = 0; rank < N_proc; rank++ ) {
        db->putScalar( "commRank", rank, AMP::Units(), AMP::Database::Check::Overwrite );
        // Create the mesh
        auto mesh0 = AMP::Mesh::Mesh::buildMesh( params );
        auto mesh  = std::dynamic_pointer_cast<AMP::Mesh::PureLogicalMesh>( mesh0 );
        AMP_ASSERT( mesh );
        // Convert an element box to a node box
        auto getNodeBox = [mesh]( AMP::Mesh::BoxMesh::Box box ) {
            auto global   = mesh->getGlobalBox();
            auto periodic = mesh->periodic();
            periodic.resize( static_cast<int>( mesh->getGeomType() ), true );
            for ( int d = 0; d < 3; d++ ) {
                if ( box.last[d] == global.last[d] && !periodic[d] )
                    box.last[d]++;
            }
            return box;
        };
        // Get the box and check the global size
        globalElemBox = mesh->getGlobalBox();
        globalNodeBox = getNodeBox( globalElemBox );
        AMP_ASSERT( mesh->numGlobalElements( AMP::Mesh::GeomType::Volume ) ==
                    globalElemBox.size().length() );
        AMP_ASSERT( mesh->numGlobalElements( AMP::Mesh::GeomType::Vertex ) ==
                    globalNodeBox.size().length() );
        // Check if we have any empty ranks and the minimum (used) box size
        if ( mesh->numLocalElements( AMP::Mesh::GeomType::Volume ) == 0 ) {
            printf( "Proc %i is empty for %i processors: (%i,%i,%i)\n",
                    rank + 1,
                    N_proc,
                    size[0],
                    size[1],
                    size[2] );
        } else {
            auto np     = mesh->numBlocks();
            auto local  = mesh->getLocalBox().size();
            auto global = globalElemBox.size();
            for ( int d = 0; d < 3; d++ ) {
                AMP_ASSERT( local[d] >= ( global[d] / np[d] ) );
                AMP_ASSERT( local[d] <= ( ( global[d] + np[d] - 1 ) / np[d] ) );
            }
        }
        // Check the iterators
        auto localElemBox = mesh->getLocalBox();
        auto localNodeBox = getNodeBox( mesh->getLocalBox() );
        size_t N_elem2    = localElemBox.size().length();
        size_t N_node2    = localNodeBox.size().length();
        AMP_ASSERT( mesh->numLocalElements( AMP::Mesh::GeomType::Vertex ) == N_node2 );
        AMP_ASSERT( mesh->numLocalElements( AMP::Mesh::GeomType::Volume ) == N_elem2 );
        AMP_ASSERT( mesh->getIterator( AMP::Mesh::GeomType::Vertex, 0 ).size() == N_node2 );
        AMP_ASSERT( mesh->getIterator( AMP::Mesh::GeomType::Volume, 0 ).size() == N_elem2 );
        // Loop through the global nodes to make sure we can correctly convert all node indices
        for ( int k = 0; k <= globalNodeBox.last[2]; k++ ) {
            for ( int j = 0; j <= globalNodeBox.last[1]; j++ ) {
                for ( int i = 0; i <= globalNodeBox.last[0]; i++ ) {
                    AMP::Mesh::BoxMesh::MeshElementIndex id1(
                        AMP::Mesh::GeomType::Vertex, 0, i, j, k );
                    auto id2 = mesh->convert( id1 );
                    auto id3 = mesh->convert( id2 );
                    AMP_ASSERT( id1 == id3 );
                    AMP_ASSERT( id2.type() == AMP::Mesh::GeomType::Vertex );
                    if ( id2.is_local() )
                        N_node++;
                    bool local = i >= localNodeBox.first[0] && i <= localNodeBox.last[0] &&
                                 j >= localNodeBox.first[1] && j <= localNodeBox.last[1] &&
                                 k >= localNodeBox.first[2] && k <= localNodeBox.last[2];
                    if ( id2.is_local() != local ) {
                        std::cout << id1 << std::endl;
                        std::cout << id2 << std::endl;
                        std::cout << id3 << std::endl;
                    }
                    AMP_ASSERT( id2.is_local() == local );
                }
            }
        }
        // Loop through the global elements to make sure we can correctly convert all node indices
        for ( int k = 0; k <= globalElemBox.last[2]; k++ ) {
            for ( int j = 0; j <= globalElemBox.last[1]; j++ ) {
                for ( int i = 0; i <= globalElemBox.last[0]; i++ ) {
                    AMP::Mesh::BoxMesh::MeshElementIndex id1(
                        AMP::Mesh::GeomType::Volume, 0, i, j, k );
                    auto id2 = mesh->convert( id1 );
                    auto id3 = mesh->convert( id2 );
                    AMP_ASSERT( id1 == id3 );
                    AMP_ASSERT( id2.type() == AMP::Mesh::GeomType::Volume );
                    if ( id2.is_local() )
                        N_elem++;
                    bool local = i >= localElemBox.first[0] && i <= localElemBox.last[0] &&
                                 j >= localElemBox.first[1] && j <= localElemBox.last[1] &&
                                 k >= localElemBox.first[2] && k <= localElemBox.last[2];
                    AMP_ASSERT( id2.is_local() == local );
                }
            }
        }
    }
    AMP_ASSERT( N_elem == globalElemBox.size().length() );
    AMP_ASSERT( N_node == globalNodeBox.size().length() );
    return true;
}


// Main function
int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );

    // Run conversion tests
    bool pass = testConversion();

    // Run different load balance sizes
    std::array<int, 3> size[] = { { 88, 1, 1 }, { 128, 1, 1 }, { 2500, 1, 1 }, { 7, 13, 17 } };
    int N_procs[]             = { 2, 3, 4, 5, 8, 12, 13, 16, 64, 87 };
    for ( auto s : size ) {
        for ( auto N : N_procs ) {
            pass = pass && run( s, { 0, 0, 0 }, N );
        }
    }

    // Run some specific problems
    pass = pass && run( { 20, 20, 20 }, { 1, 0, 0 }, 2 );

    // Shutdown AMP
    if ( pass )
        std::cout << "All tests passed\n";
    AMP::AMPManager::shutdown();
    return pass ? 0 : 1;
}
