#include <sstream>
#include <string>

#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"

#include "AMP/ampmesh/Mesh.h"
#include "AMP/utils/Writer.h"


void createMesh( AMP::UnitTest *ut, const std::string &input_file )
{
    // Create a global communicator
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Read the input file
    auto input_db = AMP::Database::parseInputFile( input_file );

    // Get the Mesh database and create the mesh parameters
    std::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( database ) );
    params->setComm( globalComm );

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh( params );

    // Create the silo writer and write the mesh
    AMP::Utilities::recursiveMkdir( "output" );
    std::string output_file                       = "output/" + input_file;
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    if ( siloWriter != nullptr ) {
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
