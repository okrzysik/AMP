#include <string>

#include "AMP/ampmesh/Mesh.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"


int main( int argc, char **argv )
{
    // Startup
    AMP::AMPManager::startup( argc, argv );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Loop through the input arguments
    for ( int i = 1; i < argc; i++ ) {

        // Get the input filename
        std::string filename = argv[i];
        auto filename2       = filename.substr( 0, filename.rfind( '.' ) );
        if ( filename2.rfind( '/' ) != std::string::npos )
            filename2 = filename2.substr( filename2.rfind( '/' ) + 1 );
        if ( filename2.rfind( '\\' ) != std::string::npos )
            filename2 = filename2.substr( filename2.rfind( '\\' ) + 1 );
        std::cout << filename2 << ": " << filename << std::endl;

        // Create the input mesh database
        auto db = std::make_shared<AMP::Database>( "Mesh" );
        db->putScalar( "MeshType", "AMP" );
        db->putScalar( "MeshName", filename2 );
        db->putScalar( "FileName", filename );

        // Create the mesh
        auto params = std::make_shared<AMP::Mesh::MeshParameters>( db );
        params->setComm( globalComm );
        auto mesh = AMP::Mesh::Mesh::buildMesh( params );

        // Create the silo writer and register the data
        auto siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
        siloWriter->registerMesh( mesh, 1 );
        globalComm.barrier();

        // Write the output
        siloWriter->setDecomposition( 1 );
        siloWriter->writeFile( filename2, 0 );
        globalComm.barrier();
    }

    // Shutdown
    AMP::AMPManager::shutdown();
    return 0;
}
