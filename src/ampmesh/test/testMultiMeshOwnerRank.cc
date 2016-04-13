#include <ampmesh/Mesh.h>
#include <utils/AMPManager.h>
#include <utils/AMP_MPI.h>
#include <utils/InputDatabase.h>
#include <utils/InputManager.h>
#include <utils/PIO.h>
#include <utils/UnitTest.h>

void testMultiMeshOwnerRank( AMP::UnitTest& ut )
{
//    std::string const exeName = "testDTKConstruction";
    std::string const exeName = "testMultiMeshOwnerRank";
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
    if ( arrayMesh )
        it = arrayMesh->getBoundaryIDIterator(AMP::Mesh::Vertex,0);
    auto arrayBoundaryMesh = mesh->Subset(it);

    // Iterate through the vertices on the boundaries of array and see if the
    // ranks are correct.
    bool failure = false;
    if ( arrayBoundaryMesh )
    {
	int bnd_comm_rank = arrayBoundaryMesh->getComm().getRank();
	bool owner_rank_is_comm_rank = true;
	it = arrayBoundaryMesh->getIterator( AMP::Mesh::Vertex, 0 );
	auto it_begin = it.begin();
	auto it_end = it.end();
	for ( it = it_begin; it != it_end; ++it )
	{
	    // If the owner rank of the vertex is the same as the comm rank of
	    // the boundary mesh that we got the vertex from then it should be
	    // locally owned.
	    owner_rank_is_comm_rank =
		( bnd_comm_rank == it->globalID().owner_rank() );

	    // If the vertex thinks it is locally owned but its owner rank
	    // and mesh comm rank dont match then this is a failure.
	    failure = (owner_rank_is_comm_rank !=
			   it->globalID().is_local() );

	    // Exit the for loop on failure.
	    if ( failure ) break;
	}
    }

    // Return pass/fail.
    if ( !failure )
    {
	ut.passes("Owner ranks are correct");
    }
    else
    {
	ut.failure( "Owner ranks failed" );
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
