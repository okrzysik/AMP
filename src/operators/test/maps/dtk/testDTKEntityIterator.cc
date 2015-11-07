
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>
#include <cstdlib>

#include "utils/shared_ptr.h"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"

#include "operators/map/dtk/DTKAMPMeshEntityIterator.h"

bool selectAll( DataTransferKit::Entity entity )
{
    return true;
}

void myTest(AMP::UnitTest *ut)
{
    std::string exeName("testDTKEntityIterator");
    std::string log_file = "output_" + exeName;
    std::string msgPrefix;
    AMP::PIO::logOnlyNodeZero(log_file);

    AMP::pout<<"Loading the  mesh"<<std::endl;
    AMP::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

    std::string input_file = "input_" + exeName;
    AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
    input_db->printClassData(AMP::plog);

    AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
    AMP::shared_ptr<AMP::Database> meshDatabase = input_db->getDatabase("Mesh");

    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(meshDatabase));
    meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh(meshParams);

    AMP::Mesh::MeshIterator mesh_iterator = mesh->getIterator( AMP::Mesh::Volume );

    DataTransferKit::EntityIterator dtk_iterator =
	AMP::Operator::AMPMeshEntityIterator( mesh_iterator, selectAll );

    // Check the size of the iterator.
    AMP_ASSERT( dtk_iterator.size() == mesh_iterator.size() );

    for ( dtk_iterator = dtk_iterator.begin(),
    	 mesh_iterator = mesh_iterator.begin();
    	  dtk_iterator != dtk_iterator.end();
    	  ++dtk_iterator, ++mesh_iterator )
    {
    	// Check the id.
    	unsigned int tmp = 0x00000000;
    	int owner_rank = mesh_iterator->globalID().owner_rank();
    	tmp += (0x007FFFFF&owner_rank) << 8;
    	AMP::Mesh::GeomType type = mesh_iterator->globalID().type();
    	tmp += ((unsigned char) type);
    	unsigned int local_id = mesh_iterator->globalID().local_id();
    	DataTransferKit::EntityId element_id = 
    	    (((AMP::Mesh::uint64)tmp)<<32) + ((AMP::Mesh::uint64)local_id);
    	AMP_ASSERT( dtk_iterator->id() == element_id );

    	// Check the entity.
    	AMP_ASSERT( dtk_iterator->entityType() == DataTransferKit::ENTITY_TYPE_VOLUME );
    	AMP_ASSERT( (unsigned) dtk_iterator->ownerRank() ==
    		    mesh_iterator->globalID().owner_rank() );
    	AMP_ASSERT( dtk_iterator->physicalDimension() == 3 );
	
    	// Check the bounding box.
    	std::vector<AMP::Mesh::MeshElement> vertices =
    	    mesh_iterator->getElements( AMP::Mesh::Vertex );
    	AMP_ASSERT( 8 == vertices.size() );
    	std::vector<double> box(6);
    	Teuchos::Tuple<double,6> element_box;
    	dtk_iterator->boundingBox( element_box );
    	for ( int i = 0; i < 8; ++i )
    	{
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
    	for ( unsigned i = 0; i < block_ids.size(); ++i )
    	{
    	    AMP_ASSERT( dtk_iterator->inBlock(block_ids[i]) 
    			== mesh_iterator->isInBlock(block_ids[i]) );
    	}
    	std::vector<int> boundary_ids = mesh->getBoundaryIDs();
    	for ( unsigned i = 0; i < boundary_ids.size(); ++i )
    	{
    	    AMP_ASSERT( dtk_iterator->onBoundary(boundary_ids[i]) ==
    			mesh_iterator->isOnBoundary(boundary_ids[i]) );
    	}
    }

    ut->passes( exeName );

}


int main(int argc, char *argv[])
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);
    AMP::UnitTest ut;

    myTest(&ut);

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   


