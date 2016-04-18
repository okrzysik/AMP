
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <map>

#include "utils/shared_ptr.h"

#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"

#include "operators/map/dtk/DTKAMPMeshEntity.h"

void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testDTKEntity" );
    std::string log_file = "output_" + exeName;
    std::string msgPrefix;
    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::pout << "Loading the  mesh" << std::endl;
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    std::string input_file = "input_" + exeName;
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> meshDatabase = input_db->getDatabase( "Mesh" );

    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( meshDatabase ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh( meshParams );

    // get iterators
    AMP::Mesh::MeshIterator vol_iterator = mesh->getIterator( AMP::Mesh::Volume );

    // map the volume ids to dtk ids
    int counter = 0;
    std::map<AMP::Mesh::MeshElementID,DataTransferKit::EntityId> vol_id_map;
    for ( vol_iterator = vol_iterator.begin();
	  vol_iterator != vol_iterator.end();
          ++vol_iterator, ++counter )
    {
	vol_id_map.emplace( vol_iterator->globalID(), counter );
    }
    int comm_rank = globalComm.getRank();    
    int comm_size = globalComm.getSize();
    std::vector<std::size_t> offsets( comm_size, 0 );
    globalComm.allGather( vol_id_map.size(), offsets.data() );
    for ( int n = 1; n < comm_size; ++n )
    {
	offsets[n] += offsets[n-1];
    }
    if ( comm_rank > 0 )
    {
	for ( auto& i : vol_id_map ) i.second += offsets[comm_rank-1];
    }

    // build the rank map.
    std::unordered_map<int,int> rank_map;
    auto global_ranks = mesh->getComm().globalRanks();
    for ( int n = 0; n < comm_size; ++n )
    {
	rank_map.emplace( global_ranks[n], n );
    }

    for ( vol_iterator = vol_iterator.begin(); vol_iterator != vol_iterator.end();
          ++vol_iterator ) {
        // Create a dtk entity.
        DataTransferKit::Entity dtk_entity =
	    AMP::Operator::AMPMeshEntity( *vol_iterator, rank_map, vol_id_map );

        // Check the id.
        uint32_t mesh_id = vol_iterator->globalID().meshID().getData();
        unsigned int owner_rank = vol_iterator->globalID().owner_rank();
        unsigned int local_id = vol_iterator->globalID().local_id();
        DataTransferKit::EntityId element_id =
	    vol_id_map.find( vol_iterator->globalID() )->second;
        AMP_ASSERT( dtk_entity.id() == element_id );

        // Check the entity.
        AMP_ASSERT( dtk_entity.topologicalDimension() == 3 );
        AMP_ASSERT( (unsigned) dtk_entity.ownerRank() == vol_iterator->globalID().owner_rank() );
        AMP_ASSERT( dtk_entity.physicalDimension() == 3 );

        // Check the bounding box.
        std::vector<AMP::Mesh::MeshElement> vertices =
            vol_iterator->getElements( AMP::Mesh::Vertex );
        AMP_ASSERT( 8 == vertices.size() );
        std::vector<double> box( 6 );
        Teuchos::Tuple<double, 6> element_box;
        dtk_entity.boundingBox( element_box );
        for ( int i = 0; i < 8; ++i ) {
            std::vector<double> coords = vertices[i].coord();
            AMP_ASSERT( coords.size() == 3 );
            AMP_ASSERT( coords[0] >= element_box[0] );
            AMP_ASSERT( coords[0] <= element_box[3] );
            AMP_ASSERT( coords[1] >= element_box[1] );
            AMP_ASSERT( coords[1] <= element_box[4] );
            AMP_ASSERT( coords[2] >= element_box[2] );
            AMP_ASSERT( coords[2] <= element_box[5] );
        }

        // Check the block/boundary.
        std::vector<int> block_ids = mesh->getBlockIDs();
        for ( unsigned i = 0; i < block_ids.size(); ++i ) {
            AMP_ASSERT( dtk_entity.inBlock( block_ids[i] ) ==
                        vol_iterator->isInBlock( block_ids[i] ) );
        }
        std::vector<int> boundary_ids = mesh->getBoundaryIDs();
        for ( unsigned i = 0; i < boundary_ids.size(); ++i ) {
            AMP_ASSERT( dtk_entity.onBoundary( boundary_ids[i] ) ==
                        vol_iterator->isOnBoundary( boundary_ids[i] ) );
        }
    }

    ut->passes( exeName );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;

    myTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
