#include <ampmesh/Mesh.h>
#include <utils/AMPManager.h>
#include <utils/AMP_MPI.h>
#include <utils/InputDatabase.h>
#include <utils/InputManager.h>
#include <utils/PIO.h>
#include <utils/UnitTest.h>

#include <algorithm>
#include <functional>
#include <unordered_map>

union id_mask {
    AMP::Mesh::uint64 id_i[2];
    char id_c[sizeof( AMP::Mesh::uint64[2] )];
};

auto hash_id( AMP::Mesh::MeshElementID id ) -> std::hash<std::string>::result_type
{
    id_mask m;

    // get the first 64 bits from the mesh id
    m.id_i[0] = id.meshID().getData();

    // construct the next 64 from the element id.
    unsigned int tmp = 0x00000000;
    if ( id.is_local() )
        tmp = 0x80000000;
    tmp += ( 0x007FFFFF & id.owner_rank() ) << 8;
    char type = (char) id.type();
    tmp += ( (unsigned char) type );
    m.id_i[1] = ( ( (AMP::Mesh::uint64) tmp ) << 32 ) + ( (AMP::Mesh::uint64) id.local_id() );

    // hash the id
    std::hash<std::string> hasher;
    return hasher( std::string( m.id_c, sizeof( AMP::Mesh::uint64[2] ) ) );
}

void makeCommRankMap( const AMP::AMP_MPI &comm, std::unordered_map<int, int> &rank_map )
{
    auto global_ranks = comm.globalRanks();
    int size          = comm.getSize();
    for ( int n = 0; n < size; ++n ) {
        rank_map.emplace( global_ranks[n], n );
    }
}

int mapElementOwnerRank( const std::unordered_map<int, int> &rank_map,
                         AMP::Mesh::MeshElement element )
{
    return rank_map.find( element.globalOwnerRank() )->second;
}

void testMultiMeshOwnerRank( AMP::UnitTest &ut )
{
    std::string const exeName   = "testMultiMeshOwnerRank";
    std::string const inputFile = "input_" + exeName;
    std::string const logFile   = "output_" + exeName;

    AMP::PIO::logAllNodes( logFile );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Parse the input File
    AMP::shared_ptr<AMP::InputDatabase> inputDatabase( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( inputFile, inputDatabase );

    // Read the mesh database
    AMP::Database::shared_ptr meshDatabase = inputDatabase->getDatabase( "Mesh" );

    // Build the mesh
    AMP::Mesh::MeshParameters::shared_ptr meshParams(
        new AMP::Mesh::MeshParameters( meshDatabase ) );
    meshParams->setComm( globalComm );
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh( meshParams );

    // Subset the mesh boundary on surface 0
    AMP::Mesh::Mesh::shared_ptr arrayMesh = mesh->Subset( "Mesh1" );
    AMP::Mesh::MeshIterator it;
    AMP::Mesh::Mesh::shared_ptr arrayBoundaryMesh;
    if ( arrayMesh ) {
        it                = arrayMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 0 );
        arrayBoundaryMesh = arrayMesh->Subset( it );
    }

    bool failure = false;
    if ( arrayBoundaryMesh ) {
        // Create the rank mapping.
        std::unordered_map<int, int> rank_map;
        makeCommRankMap( arrayBoundaryMesh->getComm(), rank_map );

        // Iterate through the vertices on the boundaries of array and see if the
        // ranks are correct.
        int bnd_comm_rank            = arrayBoundaryMesh->getComm().getRank();
        bool owner_rank_is_comm_rank = true;
        it                           = arrayBoundaryMesh->getIterator( AMP::Mesh::Vertex, 0 );
        auto it_begin                = it.begin();
        auto it_end                  = it.end();
        for ( it = it_begin; it != it_end; ++it ) {
            // If the owner rank of the vertex is the same as the comm rank of
            // the boundary mesh that we got the vertex from then it should be
            // locally owned.
            owner_rank_is_comm_rank = ( bnd_comm_rank == mapElementOwnerRank( rank_map, *it ) );

            // If the vertex thinks it is locally owned but its owner rank
            // and mesh comm rank dont match then this is a failure.
            failure = ( owner_rank_is_comm_rank != it->globalID().is_local() );

            // Exit the for loop on failure.
            if ( failure )
                break;
        }
    }

    // Return pass/fail.
    if ( !failure ) {
        ut.passes( "Owner ranks are correct" );
    } else {
        ut.failure( "Owner ranks failed" );
    }

    // Do a reduction to make sure we only get one instance of locally owned elements.
    std::vector<unsigned long long> local_ids( 0 );
    if ( arrayBoundaryMesh ) {
        it              = arrayMesh->getBoundaryIDIterator( AMP::Mesh::Volume, 0 );
        auto volBndMesh = arrayMesh->Subset( it );
        it              = volBndMesh->getIterator( AMP::Mesh::Volume, 0 );
        auto it_begin   = it.begin();
        auto it_end     = it.end();
        for ( it = it_begin; it != it_end; ++it ) {
            local_ids.push_back( hash_id( it->globalID() ) );
        }
    }
    std::vector<int> num_ids( globalComm.getSize(), -1 );
    globalComm.allGather( static_cast<int>( local_ids.size() ), num_ids.data() );
    int num_global = 0;
    for ( auto i : num_ids )
        num_global += i;
    std::vector<unsigned long long> global_ids( num_global );
    globalComm.allGather( local_ids.data(), local_ids.size(), global_ids.data() );
    failure = false;
    for ( auto i : global_ids ) {
        auto count = std::count( global_ids.begin(), global_ids.end(), i );
        if ( count > 1 ) {
            failure = true;
            break;
        }
    }

    // Return pass/fail.
    if ( !failure ) {
        ut.passes( "Global IDs are correct" );
    } else {
        ut.failure( "Repeated global ids" );
    }
}

// Main function
int main( int argc, char **argv )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;

    // Run the MultiMesh subest test
    testMultiMeshOwnerRank( ut );

    // Print the results and return
    ut.report();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
