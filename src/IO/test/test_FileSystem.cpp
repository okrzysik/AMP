#include "AMP/IO/FileSystem.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"

#include <thread>


// Helper function to record pass/failure
#define PASS_FAIL( test, MSG ) \
    do {                       \
        if ( test )            \
            ut.passes( MSG );  \
        else                   \
            ut.failure( MSG ); \
    } while ( 0 )


/********************************************************
 *  Test filesystem                                      *
 ********************************************************/
int main( int argc, char *argv[] )
{
    AMP::AMP_MPI::start_MPI( argc, argv );
    AMP::UnitTest ut;

    // Test deleting and checking if a file exists
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    if ( globalComm.getRank() == 0 ) {
        FILE *fid = fopen( "testDeleteFile.txt", "w" );
        fputs( "Temporary test", fid );
        fclose( fid );
        PASS_FAIL( AMP::IO::exists( "testDeleteFile.txt" ), "File exists" );
        AMP::IO::deleteFile( "testDeleteFile.txt" );
        PASS_FAIL( !AMP::IO::exists( "testDeleteFile.txt" ), "File deleted" );
    }

    // Test creating/deleting directories
    using namespace std::chrono_literals;
    AMP::IO::recursiveMkdir( "." );
    AMP::IO::recursiveMkdir( "testUtilitiesDir/a/b" );
    globalComm.barrier();
    std::this_thread::sleep_for( 10ms );
    PASS_FAIL( AMP::IO::exists( "testUtilitiesDir/a/b" ), "Create directory" );
    globalComm.barrier();
    if ( globalComm.getRank() == 0 ) {
        AMP::IO::deleteFile( "testUtilitiesDir/a/b" );
        AMP::IO::deleteFile( "testUtilitiesDir/a" );
        AMP::IO::deleteFile( "testUtilitiesDir" );
        std::this_thread::sleep_for( 10ms );
    }
    globalComm.barrier();
    PASS_FAIL( !AMP::IO::exists( "testUtilitiesDir/a/b" ), "Destroy directory" );

    ut.report();
    int num_failed = ut.NumFailGlobal();

    AMP::AMP_MPI::stop_MPI();
    return num_failed;
}


