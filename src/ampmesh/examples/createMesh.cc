#include <sstream>
#include <string>

#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "ampmesh/Mesh.h"
#include "utils/Writer.h"


void createMesh( AMP::UnitTest *ut, const std::string &input_file )
{
    // Create a global communicator
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Read the input file
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );

    // Get the Mesh database and create the mesh parameters
    AMP::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( database ) );
    params->setComm( globalComm );

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh( params );

    // Create the silo writer and write the mesh
    AMP::Utilities::recursiveMkdir( "output" );
    std::string output_file                       = "output/" + input_file;
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    if ( siloWriter != NULL ) {
        siloWriter->registerMesh( mesh, 1 );
        siloWriter->setDecomposition( 1 );
        siloWriter->writeFile( output_file, 0 );
    }

    // Finished
    ut->passes( "example ran to completion" );
}


int main( int argc, char **argv )
{
    // Initialize AMP
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    // Run the example
    AMP_INSIST( argc == 2, "Usage is:  createmesh  inputFile" );
    std::string filename( argv[1] );
    createMesh( &ut, filename );

    // Report the results and exit
    ut.report();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
