#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"


void testMultiMeshSubset( AMP::UnitTest &ut )
{
    std::string const exeName   = "testMultiMeshSubset";
    std::string const inputFile = "input_" + exeName;
    std::string const logFile   = "output_" + exeName;

    AMP::PIO::logAllNodes( logFile );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Parse the input File
    auto inputDatabase = AMP::Database::parseInputFile( inputFile );

    // Read the mesh database
    auto meshDatabase = inputDatabase->getDatabase( "Mesh" );

    // Build the mesh
    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( meshDatabase );
    meshParams->setComm( globalComm );
    auto mesh = AMP::Mesh::Mesh::buildMesh( meshParams );

    // Subset the mesh
    auto fooMesh = mesh->Subset( "Foo" );
    AMP::Mesh::MeshIterator it;
    if ( fooMesh != nullptr )
        it = fooMesh->getBoundaryIDIterator( AMP::Mesh::GeomType::Volume, 0 );
    auto fooBoundaryMesh = mesh->Subset( it );

    // Check the number of elements in the subset
    size_t N_local = 0;
    if ( fooBoundaryMesh != nullptr )
        N_local = fooBoundaryMesh->numLocalElements( AMP::Mesh::GeomType::Volume );
    size_t N_global = globalComm.sumReduce( N_local );
    if ( N_global == 12 )
        ut.passes( "Subset worked correctly" );
    else
        ut.failure( "Subset failed" );
}

// Main function
int main( int argc, char **argv )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;

    // Run the MultiMesh subest test
    testMultiMeshSubset( ut );

    // Print the results and return
    ut.report();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
