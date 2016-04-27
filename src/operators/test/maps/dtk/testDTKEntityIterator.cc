
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

#include "operators/map/dtk/DTKAMPMeshEntityIterator.h"

bool selectAll( DataTransferKit::Entity entity ) { return true; }

void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testDTKEntityIterator" );
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

    // get the volume iterator
    AMP::Mesh::MeshIterator vol_iterator = mesh->getIterator( AMP::Mesh::Volume );

    // map the volume ids to dtk ids
    int counter = 0;
    AMP::shared_ptr<std::map<AMP::Mesh::MeshElementID,DataTransferKit::EntityId> >vol_id_map
	= std::make_shared<std::map<AMP::Mesh::MeshElementID,DataTransferKit::EntityId> >();
    for ( vol_iterator = vol_iterator.begin();
	  vol_iterator != vol_iterator.end();
          ++vol_iterator, ++counter )
    {
	vol_id_map->emplace( vol_iterator->globalID(), counter );
    }
    int comm_rank = globalComm.getRank();    
    int comm_size = globalComm.getSize();
    std::vector<std::size_t> offsets( comm_size, 0 );
    globalComm.allGather( vol_id_map->size(), offsets.data() );
    for ( int n = 1; n < comm_size; ++n )
    {
	offsets[n] += offsets[n-1];
    }
    if ( comm_rank > 0 )
    {
	for ( auto& i : *vol_id_map ) i.second += offsets[comm_rank-1];
    }
    
    // make the rank map.
    AMP::shared_ptr<std::unordered_map<int,int> > rank_map =
	std::make_shared<std::unordered_map<int,int> >();    
    auto global_ranks = mesh->getComm().globalRanks();
    int size = mesh->getComm().getSize();
    for ( int n = 0; n < size; ++n )
    {
	rank_map->emplace( global_ranks[n], n );
    }
    
    DataTransferKit::EntityIterator dtk_iterator =
        AMP::Operator::AMPMeshEntityIterator( rank_map, vol_id_map, vol_iterator, selectAll );

    // Check the size of the iterator.
    AMP_ASSERT( dtk_iterator.size() == vol_iterator.size() );

    for ( dtk_iterator = dtk_iterator.begin(), vol_iterator = vol_iterator.begin();
          dtk_iterator != dtk_iterator.end();
          ++dtk_iterator, ++vol_iterator ) {
        // Check the id.
        DataTransferKit::EntityId element_id =
	    vol_id_map->find( vol_iterator->globalID() )->second;	    
        AMP_ASSERT( dtk_iterator->id() == element_id );

        // Check the entity.
        AMP_ASSERT( dtk_iterator->topologicalDimension() == 3 );
        AMP_ASSERT( (unsigned) dtk_iterator->ownerRank() ==
                    vol_iterator->globalID().owner_rank() );
        AMP_ASSERT( dtk_iterator->physicalDimension() == 3 );

        // Check the bounding box.
        std::vector<AMP::Mesh::MeshElement> vertices =
            vol_iterator->getElements( AMP::Mesh::Vertex );
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
