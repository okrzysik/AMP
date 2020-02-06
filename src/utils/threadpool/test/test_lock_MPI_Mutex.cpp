#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/threadpool/Mutex.h"
#include "AMP/utils/threadpool/ThreadPool.h"

#include "ProfilerApp.h"

#include <vector>


volatile int global_count  = 0;
volatile bool global_start = false;
AMP::Mutex global_lock( true );


void test_lock( AMP::AMP_MPI comm, int N, bool call_sleep )
{
    while ( !global_start )
        sched_yield();
    for ( int i = 0; i < N; i++ ) {
        // Acquire the lock
        AMP::lock_MPI_Mutex( global_lock, comm );
        PROFILE_START( "work", 2 );
        comm.barrier();
        // Check and increment count
        int tmp = global_count++;
        if ( tmp != 0 )
            AMP_ERROR( "Invalid count" );
        // Acquire the lock a second time, then release
        global_lock.lock();
        global_lock.unlock();
        // Sleep for a while
        sched_yield();
        if ( call_sleep )
            AMP::Utilities::sleep_ms( 20 );
        // Check and decrement count
        tmp = global_count--;
        if ( tmp != 1 )
            AMP_ERROR( "Invalid count" );
        // Release the mutex
        PROFILE_STOP( "work", 2 );
        global_lock.unlock();
        // Try to add some random waits
        for ( int j = 0; j < rand() % 10; j++ ) {
            sched_yield();
            timespec duration;
            duration.tv_sec  = 0;
            duration.tv_nsec = 100000 * ( rand() % 5 );
            nanosleep( &duration, nullptr );
        }
    }
}


int main( int argc, char *argv[] )
{
    // Initialize AMP
    AMP::AMPManager::startup( argc, argv );
    PROFILE_ENABLE( 2 );
    PROFILE_ENABLE_TRACE();
    PROFILE_START( "main" );

    {
        // Create the thread pool
        int N_threads = 8;
        AMP::ThreadPool tpool( N_threads );
        AMP::AMP_MPI comm_world( MPI_COMM_WORLD );
        comm_world.barrier();

        // Check the duration of the sleep functions
        double t0 = AMP::AMP_MPI::time();
        AMP::Utilities::sleep_ms( 1 );
        double t1 = AMP::AMP_MPI::time();
        std::cout << "AMP::Utilities::sleep_ms(1) = " << t1 - t0 << std::endl;

        // Run a single lock test
        AMP::pout << "Running single lock test\n";
        PROFILE_START( "single" );
        global_start = false;
        std::vector<AMP::ThreadPoolID> ids;
        for ( int i = 0; i < N_threads; i++ )
            ids.push_back( TPOOL_ADD_WORK( &tpool, test_lock, ( comm_world.dup(), 1, true ) ) );
        global_start = true;
        tpool.wait_all( ids );
        ids.clear();
        comm_world.barrier();
        PROFILE_STOP( "single" );

        // Run multiple lock tests
        AMP::pout << "Running multiple lock test\n";
        PROFILE_START( "multiple" );
        global_start = false;
        int N_it     = 100;
        double start = AMP::AMP_MPI::time();
        for ( int i = 0; i < N_threads; i++ )
            ids.push_back( TPOOL_ADD_WORK( &tpool, test_lock, ( comm_world.dup(), N_it, false ) ) );
        global_start = true;
        tpool.wait_all( ids );
        ids.clear();
        comm_world.barrier();
        double stop = AMP::AMP_MPI::time();
        PROFILE_STOP( "multiple" );
        AMP::pout << "   Time to acquire global MPI lock was " << ( stop - start ) / N_it
                  << " seconds/iteration\n";
    }

    // Finalize
    AMP::pout << "Test ran sucessfully\n";
    PROFILE_STOP( "main" );
    PROFILE_SAVE( "test_lock_MPI_Mutex" );
    AMP::AMPManager::shutdown();
    return 0;
}
