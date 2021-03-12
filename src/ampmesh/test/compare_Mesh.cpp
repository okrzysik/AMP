#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/ampmesh/structured/BoxMesh.h"
#include "AMP/ampmesh/testHelpers/meshTests.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Array.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"


// Compare meshes
bool compare( const std::string &filename )
{
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Create the meshes
    auto input_db = AMP::Database::parseInputFile( filename );
    std::vector<std::shared_ptr<AMP::Mesh::Mesh>> meshes;
    for ( int i = 1; i < 1000; i++ ) {
        auto key = "Mesh_" + std::to_string( i );
        if ( !input_db->keyExists( key ) )
            break;
        auto db     = input_db->getDatabase( key );
        auto params = std::make_shared<AMP::Mesh::MeshParameters>( db );
        params->setComm( globalComm );
        meshes.push_back( AMP::Mesh::Mesh::buildMesh( params ) );
    }

    // Compare the meshes
    AMP::Array<int> result( meshes.size(), meshes.size() );
    result.fill( 0 );
    for ( size_t i = 0; i < meshes.size(); i++ ) {
        for ( size_t j = 0; j < meshes.size(); j++ )
            result( i, j ) = AMP::Mesh::Mesh::compare( *meshes[i], *meshes[j] );
    }

    // Print the results
    if ( globalComm.getRank() == 0 ) {
        result.print( std::cout, "compare" );
        std::cout << std::endl;
    }

    // Check we get the same result on all ranks
    int error = 0;
    auto tmp  = result;
    globalComm.bcast( tmp.data(), tmp.length(), 0 );
    if ( tmp != result ) {
        std::cerr << "Results do not match across ranks:\n";
        std::cerr << "   rank:\n";
        result.print( std::cerr, "tmp", "   " );
        error = 1;
    }

    // Check the answer (if present)
    if ( !input_db->keyExists( "ans" ) ) {
        auto ans  = input_db->getArray<int>( "ans" );
        bool test = ans == result;
        if ( !test ) {
            std::cerr << "Answer does not match\n";
            error = 1;
        }
    }

    return globalComm.maxReduce( error );
}


// Main function
int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );

    if ( argc != 2 ) {
        std::cerr << "compare_Mesh input_file\n";
        AMP::AMPManager::shutdown();
        return 1;
    }

    int result = compare( argv[1] );

    AMP::AMPManager::shutdown();
    return result;
}
