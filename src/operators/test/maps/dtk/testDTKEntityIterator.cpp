#include "AMP/IO/PIO.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/operators/map/dtk/DTKAMPMeshEntityIterator.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include <memory>

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>

bool selectAll( DataTransferKit::Entity entity ) { return true; }

static void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testDTKEntityIterator" );
    std::string log_file = "output_" + exeName;
    std::string msgPrefix;
    AMP::logOnlyNodeZero( log_file );

    AMP::pout << "Loading the  mesh" << std::endl;

    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    std::string input_file = "input_" + exeName;
    auto input_db          = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto meshDatabase = input_db->getDatabase( "Mesh" );

    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( meshDatabase );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto mesh = AMP::Mesh::MeshFactory::create( meshParams );

    // get the volume iterator
    auto vol_iterator = mesh->getIterator( AMP::Mesh::GeomType::Volume );

    // map the volume ids to dtk ids
    int counter = 0;
    auto vol_id_map =
        std::make_shared<std::map<AMP::Mesh::MeshElementID, DataTransferKit::EntityId>>();
    for ( vol_iterator = vol_iterator.begin(); vol_iterator != vol_iterator.end();
          ++vol_iterator, ++counter ) {
        vol_id_map->emplace( vol_iterator->globalID(), counter );
    }
    int comm_rank = globalComm.getRank();
    int comm_size = globalComm.getSize();
    std::vector<std::size_t> offsets( comm_size, 0 );
    globalComm.allGather( vol_id_map->size(), offsets.data() );
    for ( int n = 1; n < comm_size; ++n ) {
        offsets[n] += offsets[n - 1];
    }
    if ( comm_rank > 0 ) {
        for ( auto &i : *vol_id_map )
            i.second += offsets[comm_rank - 1];
    }

    // make the rank map.
    auto rank_map     = std::make_shared<std::unordered_map<int, int>>();
    auto global_ranks = mesh->getComm().globalRanks();
    int size          = mesh->getComm().getSize();
    for ( int n = 0; n < size; ++n ) {
        rank_map->emplace( global_ranks[n], n );
    }

    auto dtk_iterator =
        AMP::Operator::AMPMeshEntityIterator( rank_map, vol_id_map, vol_iterator, selectAll );

    // Check the size of the iterator.
    AMP_ASSERT( dtk_iterator.size() == vol_iterator.size() );

    for ( dtk_iterator = dtk_iterator.begin(), vol_iterator = vol_iterator.begin();
          dtk_iterator != dtk_iterator.end();
          ++dtk_iterator, ++vol_iterator ) {
        // Check the id.
        DataTransferKit::EntityId element_id = vol_id_map->find( vol_iterator->globalID() )->second;
        AMP_ASSERT( dtk_iterator->id() == element_id );

        // Check the entity.
        AMP_ASSERT( dtk_iterator->topologicalDimension() == 3 );
        AMP_ASSERT( (unsigned) dtk_iterator->ownerRank() == vol_iterator->globalID().owner_rank() );
        AMP_ASSERT( dtk_iterator->physicalDimension() == 3 );

        // Check the bounding box.
        std::vector<AMP::Mesh::MeshElement> vertices =
            vol_iterator->getElements( AMP::Mesh::GeomType::Vertex );
        AMP_ASSERT( 8 == vertices.size() );
        std::vector<double> box( 6 );
        Teuchos::Tuple<double, 6> element_box;
        dtk_iterator->boundingBox( element_box );
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
            AMP_ASSERT( dtk_iterator->inBlock( block_ids[i] ) ==
                        vol_iterator->isInBlock( block_ids[i] ) );
        }
        std::vector<int> boundary_ids = mesh->getBoundaryIDs();
        for ( unsigned i = 0; i < boundary_ids.size(); ++i ) {
            AMP_ASSERT( dtk_iterator->onBoundary( boundary_ids[i] ) ==
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
