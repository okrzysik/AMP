#include "AMP/ampmesh/Mesh.h"
#include "AMP/operators/map/dtk/DTKAMPMeshEntitySet.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include <memory>

#include <cstdlib>
#include <iostream>
#include <string>


bool selectAll( DataTransferKit::Entity entity ) { return true; }

static void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testDTKEntitySet" );
    std::string log_file = "output_" + exeName;
    std::string msgPrefix;
    AMP::logOnlyNodeZero( log_file );

    AMP::pout << "Loading the  mesh" << std::endl;

    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    std::string input_file = "input_" + exeName;
    auto input_db          = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    std::shared_ptr<AMP::Database> meshDatabase = input_db->getDatabase( "Mesh" );

    std::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( meshDatabase ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh( meshParams );

    int gcw                               = 0;
    AMP::Mesh::MeshIterator mesh_iterator = mesh->getIterator( AMP::Mesh::GeomType::Volume, gcw );

    // Make an entity set.
    std::shared_ptr<DataTransferKit::EntitySet> dtk_entity_set(
        new AMP::Operator::AMPMeshEntitySet( mesh ) );
    AMP_ASSERT( 3 == dtk_entity_set->physicalDimension() );
    AMP_ASSERT( mesh->getComm().getRank() == dtk_entity_set->communicator()->getRank() );
    AMP_ASSERT( mesh->getComm().getSize() == dtk_entity_set->communicator()->getSize() );

    // Check the mesh with an iterator.
    DataTransferKit::EntityIterator dtk_iterator = dtk_entity_set->entityIterator( 3 );
    bool caught_exception                        = false;
    AMP_ASSERT( dtk_iterator.size() == mesh_iterator.size() );
    for ( dtk_iterator = dtk_iterator.begin(), mesh_iterator = mesh_iterator.begin();
          dtk_iterator != dtk_iterator.end();
          ++dtk_iterator, ++mesh_iterator ) {
        // Check with the iterator.
        {
            // Check the entity.
            AMP_ASSERT( dtk_iterator->topologicalDimension() == 3 );
            AMP_ASSERT( (unsigned) dtk_iterator->ownerRank() ==
                        mesh_iterator->globalID().owner_rank() );
            AMP_ASSERT( dtk_iterator->physicalDimension() == 3 );

            // Check the bounding box.
            std::vector<AMP::Mesh::MeshElement> vertices =
                mesh_iterator->getElements( AMP::Mesh::GeomType::Vertex );
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
                            mesh_iterator->isInBlock( block_ids[i] ) );
            }
            std::vector<int> boundary_ids = mesh->getBoundaryIDs();
            for ( unsigned i = 0; i < boundary_ids.size(); ++i ) {
                AMP_ASSERT( dtk_iterator->onBoundary( boundary_ids[i] ) ==
                            mesh_iterator->isOnBoundary( boundary_ids[i] ) );
            }

            // Check that we get 8 nodes from the adjacency function.
            try {
                Teuchos::Array<DataTransferKit::Entity> nodes;
                dtk_entity_set->getAdjacentEntities( *dtk_iterator, 0, nodes );
                AMP_ASSERT( 8 == nodes.size() );
            } catch ( std::runtime_error &e ) {
                caught_exception = true;
            }
        }
    }

    if ( caught_exception ) {
        ut->expected_failure( "caught getAdjacentEntities exception" );
    } else {
        ut->failure( "failed to catch getAdjacentEntities exception" );
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
