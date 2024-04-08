#include "AMP/IO/RestartManager.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/UnitTest.h"

#include "ProfilerApp.h"

#include <algorithm>
#include <complex>
#include <iostream>
#include <random>


void record( AMP::UnitTest &ut, bool pass, const std::string &msg )
{
    if ( pass )
        ut.passes( msg );
    else
        ut.failure( msg );
}


void testRestartManagerComms( AMP::UnitTest &ut )
{
    PROFILE( "MAIN" );

    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    int size = globalComm.getSize();
    int rank = globalComm.getRank();

    // Create some initial (special) comms
    std::vector<AMP::AMP_MPI> comms;
    comms.emplace_back( AMP_COMM_WORLD );
    comms.emplace_back( AMP_COMM_SELF );
    comms.emplace_back( AMP_COMM_NULL );
#ifdef AMP_USE_MPI
    comms.emplace_back( MPI_COMM_WORLD );
#else
    comms.emplace_back( AMP_COMM_WORLD );
#endif

    // Create a comm on only some ranks
    auto tmp = globalComm.split( rank < 2 ? 0 : 1 );
    if ( rank < 2 )
        comms.push_back( tmp );
    AMP::AMP_MPI tmp2( tmp.getCommunicator(), false );
    record( ut, tmp.hash() == tmp2.hash(), "hash is tied to MPI" );

    // Create some random comms
    std::mt19937 gen( globalComm.rand() );
    std::uniform_int_distribution<int> dis( 0, size );
    auto key = globalComm.globalRanks();
    for ( int i = 0; i < 20; i++ ) {
        std::shuffle( key.begin(), key.end(), gen );
        comms.push_back( globalComm.split( rank < dis( gen ) ? 0 : 1, key[rank] ) );
    }

    // Duplicate all comms
    size_t N = comms.size();
    for ( size_t i = 0; i < N; i++ )
        comms.push_back( comms[i].dup() );

    // Create the restart manager and register data
    AMP::IO::RestartManager writer;
    for ( auto &comm : comms )
        writer.registerComm( comm );

    // Write the restart data
    writer.write( "testRestartDataComm" );
    writer = AMP::IO::RestartManager();

    // Read the restart data
    AMP::IO::RestartManager reader( "testRestartDataComm" );
    std::vector<AMP::AMP_MPI> comms2;
    for ( auto &comm : comms )
        comms2.push_back( reader.getComm( comm.hash() ) );
    reader = AMP::IO::RestartManager();

    // Check the built-in comms
    record( ut, comms[0] == comms2[0], "Load AMP_COMM_WORLD" );
    record( ut, comms[1] == comms2[1], "Load AMP_COMM_SELF" );
    record( ut, comms[2] == comms2[2], "Load AMP_COMM_NULL" );
    record( ut, comms[3] == comms2[3], "Load MPI_COMM_WORLD" );

    // Check all comms are similar
    bool test = true;
    for ( size_t i = 4; i < comms.size(); i++ )
        test = test && comms[i].compare( comms2[i] ) > 0 &&
               comms[i].hashRanks() == comms2[i].hashRanks();
    record( ut, test, "Load other comms" );

    // Check all comm dups are similar
    test = true;
    for ( size_t i = 4; i < N; i++ )
        test = test && comms2[i].compare( comms2[i + N] ) > 0 && comms2[i] != comms2[i + N] &&
               comms2[i].hashRanks() == comms2[i + N].hashRanks();
    record( ut, test, "Duplicated comms are similar" );
}


int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    PROFILE_ENABLE();
    AMP::AMP_MPI::changeProfileLevel( 0 );

    testRestartManagerComms( ut );
    PROFILE_SAVE( "testRestartDataComm", 1 );

    int N_failed = ut.NumFailGlobal();
    ut.report();
    ut.reset();
    AMP::AMPManager::shutdown();
    return N_failed;
}
