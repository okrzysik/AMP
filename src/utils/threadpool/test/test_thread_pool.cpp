#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/threadpool/ThreadPool.h"

#include "ProfilerApp.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>


using namespace AMP;
using AMP::Utilities::stringf;


#define to_ns( x ) std::chrono::duration_cast<std::chrono::nanoseconds>( x ).count()
#define to_ms( x ) std::chrono::duration_cast<std::chrono::milliseconds>( x ).count()


// Function to waste CPU cycles
void waste_cpu( int N )
{
    if ( N > 10000 )
        PROFILE_START( "waste_cpu", 2 );
    double x = 1.0;
    N        = std::max( 10, N );
    {
        double pi = 3.141592653589793;
        for ( int i = 0; i < N; i++ )
            x = sqrt( x * exp( pi / x ) );
    } // style to limit gcov hits
    if ( fabs( x - 2.926064057273157 ) > 1e-12 )
        abort();
    if ( N > 10000 )
        PROFILE_STOP( "waste_cpu", 2 );
}


// Sleep for the given time
// Note: since we may encounter interrupts, we may not sleep for the desired time
//   so we need to perform the sleep in a loop
void sleep_ms( int64_t N )
{
    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::high_resolution_clock::now();
    while ( to_ms( t2 - t1 ) < N ) {
        int N2 = N - to_ms( t2 - t1 );
        std::this_thread::sleep_for( std::chrono::milliseconds( N2 ) );
        t2 = std::chrono::high_resolution_clock::now();
    }
}
void sleep_s( int N ) { sleep_ms( 1000 * N ); }


// Function to sleep for N seconds then increment a global count
static volatile int global_sleep_count = 0;
void sleep_inc( int N )
{
    PROFILE_START( "sleep_inc" );
    sleep_s( N );
    ++global_sleep_count;
    PROFILE_STOP( "sleep_inc" );
}
void sleep_inc2( double x )
{
    sleep_ms( static_cast<int>( round( x * 1000 ) ) );
    ++global_sleep_count;
}
void sleep_msg( double x, const std::string &msg )
{
    PROFILE_START( msg );
    sleep_ms( static_cast<int>( round( x * 1000 ) ) );
    NULL_USE( msg );
    PROFILE_STOP( msg );
}
bool check_inc( int N ) { return global_sleep_count == N; }


// Function to return the processor for the given thread
std::mutex print_processor_mutex;

void print_processor( ThreadPool *tpool )
{
    AMP::AMP_MPI comm( MPI_COMM_WORLD );
    int rank      = comm.getRank();
    int thread    = tpool->getThreadNumber();
    int processor = ThreadPool::getCurrentProcessor();
    auto tmp      = stringf( "%i:  Thread,proc = %i,%i\n", rank, thread, processor );
    sleep_ms( 10 * rank );
    print_processor_mutex.lock();
    pout << tmp;
    print_processor_mutex.unlock();
    sleep_ms( 100 );
}


// Function to test how a member thread interacts with the thread pool
int test_member_thread( ThreadPool *tpool )
{
    int N_errors = 0;
    // Member threads are not allowed to wait for the pool to finish
    try {
        tpool->wait_pool_finished();
        N_errors++;
    } catch ( ... ) {
    }
    // Member threads are not allowed to change the size of the pool
    try {
        tpool->wait_pool_finished();
        N_errors++;
    } catch ( ... ) {
    }
    return N_errors;
}


/******************************************************************
 * Test the TPOOL_ADD_WORK macro with variable number of arguments *
 ******************************************************************/
static int myfun0() { return 0; }
static int myfun1( int ) { return 1; }
static int myfun2( int, float ) { return 2; }
static int myfun3( int, float, double ) { return 3; }
static int myfun4( int, float, double, char ) { return 4; }
static int myfun5( int, float, double, char, std::string ) { return 5; }
static int myfun6( int, float, double, char, std::string, int ) { return 6; }
static int myfun7( int, float, double, char, std::string, int, int ) { return 7; }
static void test_function_arguments( ThreadPool *tpool, UnitTest &ut )
{
    // Test some basic types of instantiations
    bool pass = true;
    printp( "Testing arguments:\n" );
    ThreadPoolID id0 = TPOOL_ADD_WORK( tpool, myfun0, ( nullptr ) );
    ThreadPoolID id1 = TPOOL_ADD_WORK( tpool, myfun1, ( (int) 1 ) );
    ThreadPoolID id2 = TPOOL_ADD_WORK( tpool, myfun2, ( (int) 1, (float) 2 ) );
    ThreadPoolID id3 = TPOOL_ADD_WORK( tpool, myfun3, ( (int) 1, (float) 2, (double) 3 ) );
    ThreadPoolID id4 =
        TPOOL_ADD_WORK( tpool, myfun4, ( (int) 1, (float) 2, (double) 3, (char) 4 ) );
    ThreadPoolID id5 = TPOOL_ADD_WORK(
        tpool, myfun5, ( (int) 1, (float) 2, (double) 3, (char) 4, std::string( "test" ) ) );
    ThreadPoolID id52 = TPOOL_ADD_WORK(
        tpool, myfun5, ( (int) 1, (float) 2, (double) 3, (char) 4, std::string( "test" ) ), -1 );
    ThreadPoolID id6 = TPOOL_ADD_WORK(
        tpool,
        myfun6,
        ( (int) 1, (float) 2, (double) 3, (char) 4, std::string( "test" ), (int) 1 ) );
    ThreadPoolID id7 = TPOOL_ADD_WORK(
        tpool,
        myfun7,
        ( (int) 1, (float) 2, (double) 3, (char) 4, std::string( "test" ), (int) 1, (int) 1 ) );
    tpool->wait_pool_finished();
    pass = pass && tpool->isFinished( id0 );
    pass = pass && tpool->getFunctionRet<int>( id0 ) == 0;
    pass = pass && tpool->getFunctionRet<int>( id1 ) == 1;
    pass = pass && tpool->getFunctionRet<int>( id2 ) == 2;
    pass = pass && tpool->getFunctionRet<int>( id3 ) == 3;
    pass = pass && tpool->getFunctionRet<int>( id4 ) == 4;
    pass = pass && tpool->getFunctionRet<int>( id5 ) == 5;
    pass = pass && tpool->getFunctionRet<int>( id52 ) == 5;
    pass = pass && tpool->getFunctionRet<int>( id6 ) == 6;
    pass = pass && tpool->getFunctionRet<int>( id7 ) == 7;
    if ( pass )
        ut.passes( "Calling function with default arguments" );
    else
        ut.failure( "Error calling function with default arguments" );
}


/******************************************************************
 * Examples to derive a user work item                             *
 ******************************************************************/
class UserWorkItemVoid final : public ThreadPool::WorkItem
{
public:
    // User defined constructor (does not need to match any interfaces)
    explicit UserWorkItemVoid( int dummy )
    {
        // User initialized variables
        NULL_USE( dummy );
    }
    // User defined run (can do anything)
    void run() override
    {
        // Perform the tasks
        printf( "Hello work from UserWorkItem (void)" );
    }
    // Will the routine return a result
    bool has_result() const override { return false; }
    // User defined destructor
    ~UserWorkItemVoid() override = default;
};
class UserWorkItemInt final : public ThreadPool::WorkItemRet<int>
{
public:
    // User defined constructor (does not need to match any interfaces)
    explicit UserWorkItemInt( int dummy )
    {
        // User initialized variables
        NULL_USE( dummy );
    }
    // User defined run (can do anything)
    void run() override
    {
        // Perform the tasks
        printf( "Hello work from UserWorkItem (int)" );
        // Store the results (it's type will match the template)
        ThreadPool::WorkItemRet<int>::d_result = 1;
    }
    // User defined destructor
    ~UserWorkItemInt() override = default;
};


/******************************************************************
 * test the time to run N tasks in parallel                        *
 ******************************************************************/
template<class Ret, class... Args>
inline double launchAndTime( ThreadPool &tpool, int N, Ret ( *routine )( Args... ), Args... args )
{
    tpool.wait_pool_finished();
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<ThreadPoolID> ids( N );
    for ( int i = 0; i < N; i++ )
        ids[i] = ThreadPool_add_work( &tpool, 0, routine, args... );
    NULL_USE( ids );
    tpool.wait_all( ids );
    // tpool.wait_pool_finished();
    auto stop = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>( stop - start ).count();
}


// Move constructor function
volatile ThreadPoolID f1( volatile ThreadPoolID a ) { return a; }
ThreadPoolID f2( ThreadPoolID a ) { return a; }


/******************************************************************
 * Test the basic functionallity of the atomics                    *
 ******************************************************************/
template<class T>
T atomic_compare_and_swap( volatile std::atomic<T> &x, T e, T d )
{
    auto e2 = e;
    return x.compare_exchange_weak( e2, d );
}
bool test_atomics()
{
    bool pass                        = true;
    volatile std::atomic_int32_t i32 = 32;
    volatile std::atomic_int32_t i64 = 64;
    pass                             = pass && ++i32 == 33 && ++i64 == 65;
    pass                             = pass && --i32 == 32 && --i64 == 64;
    pass                             = pass && i32.fetch_add( 2 ) != 34;
    pass                             = pass && i64.fetch_add( 4 ) != 68;
    pass                             = pass && !atomic_compare_and_swap( i32, 0, 0 );
    pass                             = pass && !atomic_compare_and_swap( i64, 0, 0 );
    pass                             = pass && atomic_compare_and_swap( i32, 34, 32 );
    pass                             = pass && atomic_compare_and_swap( i64, 68, 64 );
    pass                             = pass && i32 == 32 && i64 == 64;
    return pass;
}


/******************************************************************
 * Test the move operator for thread_id                            *
 ******************************************************************/
void test_move_thread_id( UnitTest &, ThreadPool &tpool )
{
    auto id = TPOOL_ADD_WORK( &tpool, waste_cpu, ( 5 ) );
    tpool.wait( id );
    std::vector<ThreadPoolID> ids2;
    ids2.push_back( TPOOL_ADD_WORK( &tpool, waste_cpu, ( 5 ) ) );
    tpool.wait( ids2[0] );
    auto id1          = f1( id );         // move-construct from rvalue temporary
    auto id2          = std::move( id1 ); // move-construct from xvalue
    volatile auto id3 = f2( id );         // move-construct from rvalue temporary
    volatile auto id4 = std::move( id3 ); // move-construct from xvalue
    id2.reset();
    id4.reset();
}


/******************************************************************
 * Check that the threads work in parallel                         *
 ******************************************************************/
void test_sleep_parallel( UnitTest &ut, ThreadPool &tpool )
{
    int N_threads    = tpool.getNumThreads();
    int N_procs      = ThreadPool::getNumberOfProcessors();
    int N_procs_used = std::min<int>( N_procs, N_threads );
    tpool.wait_pool_finished();
    auto start = std::chrono::high_resolution_clock::now();
    sleep_inc( 1 );
    auto stop             = std::chrono::high_resolution_clock::now();
    double sleep_serial   = std::chrono::duration<double>( stop - start ).count();
    double sleep_parallel = launchAndTime( tpool, N_threads, sleep_inc, 1 );
    double sleep_speedup  = N_procs_used * sleep_serial / sleep_parallel;
    printp( "Speedup on %i sleeping threads: %0.3f\n", N_procs_used, sleep_speedup );
    printp( "   ts = %0.3f, tp = %0.3f\n", sleep_serial, sleep_parallel );
    if ( fabs( sleep_serial - 1.0 ) < 0.05 && fabs( sleep_parallel - 1.0 ) < 0.25 &&
         sleep_speedup > 3 )
        ut.passes( "Passed thread sleep" );
    else
        ut.failure( "Failed thread sleep" );
}
void test_work_parallel( UnitTest &ut, ThreadPool &tpool )
{
    int N_threads    = tpool.getNumThreads();
    int N_procs      = ThreadPool::getNumberOfProcessors();
    int N_procs_used = std::min<int>( N_procs, N_threads );
    if ( N_procs_used > 1 ) {
        AMP::AMP_MPI comm( MPI_COMM_WORLD );
        // Use a non-blocking serialization of the MPI processes
        // if we do not have a sufficient number of processors
        int rank           = comm.getRank();
        int size           = comm.getSize();
        bool serialize_mpi = N_procs < N_threads * size;
        if ( serialize_mpi && rank > 0 ) {
            int buf;
            auto request = comm.Irecv( &buf, 1, rank - 1, 0 );
            comm.wait( request );
        }
        int N = 20000000; // Enough work to keep the processor busy for ~ 1 s
        // Run in serial
        auto start = std::chrono::high_resolution_clock::now();
        waste_cpu( N );
        auto stop = std::chrono::high_resolution_clock::now();
        double ts = std::chrono::duration<double>( stop - start ).count();
        // Run in parallel
        double tp      = launchAndTime( tpool, N_procs_used, waste_cpu, N );
        double tp2     = launchAndTime( tpool, N_procs_used, waste_cpu, N / 1000 );
        double speedup = N_procs_used * ts / tp;
        printp( "Speedup on %i procs: %0.3f\n", N_procs_used, speedup );
        printp( "   ts = %0.3f, tp = %0.3f, tp2 = %0.3f\n", ts, tp, tp2 );
        if ( speedup > 1.4 ) {
            ut.passes( "Passed speedup test" );
        } else {
#ifdef USE_GCOV
            ut.expected_failure( "Times do not indicate tests are running in parallel (gcov)" );
#else
            ut.failure( "Times do not indicate tests are running in parallel" );
#endif
        }
        if ( serialize_mpi ) {
            if ( rank < size - 1 )
                comm.send( &N, 1, rank + 1, 0 );
            if ( rank == size - 1 ) {
                for ( int i = 0; i < size - 1; i++ )
                    comm.send( &N, 1, i, 1 );
            } else {
                int buf;
                auto request = comm.Irecv( &buf, 1, size - 1, 1 );
                comm.wait( request );
            }
        }
    } else {
        ut.expected_failure( "Testing thread performance with less than 1 processor" );
    }
}


/******************************************************************
 * Test adding a work item with a dependency                       *
 ******************************************************************/
void test_work_dependency( UnitTest &ut, ThreadPool &tpool )
{
    std::vector<ThreadPoolID> ids;
    ids.reserve( 5 );
    global_sleep_count = 0; // Reset the count before this test
    ThreadPoolID id0;
    auto id1    = TPOOL_ADD_WORK( &tpool, sleep_inc, ( 1 ) );
    auto id2    = TPOOL_ADD_WORK( &tpool, sleep_inc, ( 2 ) );
    auto *wait1 = ThreadPool::createWork( check_inc, 1 );
    auto *wait2 = ThreadPool::createWork( check_inc, 2 );
    wait1->add_dependency( id0 );
    wait1->add_dependency( id1 );
    wait2->add_dependency( id1 );
    wait2->add_dependency( id2 );
    ids.clear();
    ids.push_back( tpool.add_work( wait1 ) );
    ids.push_back( tpool.add_work( wait2 ) );
    tpool.wait_all( ids );
    if ( !tpool.getFunctionRet<bool>( ids[0] ) || !tpool.getFunctionRet<bool>( ids[1] ) )
        ut.failure( "Failed to wait on required dependency" );
    else
        ut.passes( "Dependencies" );
    tpool.wait_pool_finished();
    // Test waiting on more dependencies than in the thread pool (changing priorities)
    ids.clear();
    for ( size_t i = 0; i < 20; i++ )
        ids.push_back( TPOOL_ADD_WORK( &tpool, sleep_inc2, ( 0.1 ) ) );
    auto *wait3 = ThreadPool::createWork( sleep_inc2, 0 );
    wait3->add_dependencies( ids );
    auto id = tpool.add_work( wait3, 50 );
    tpool.wait( id );
    auto pass = true;
    for ( auto &tmp : ids )
        pass = pass && tmp.finished();
    ids.clear();
    if ( pass )
        ut.passes( "Dependencies2" );
    else
        ut.failure( "Dependencies2" );
    // Check that we can handle more complex dependencies
    id1 = TPOOL_ADD_WORK( &tpool, sleep_inc2, ( 0.5 ) );
    for ( int i = 0; i < 10; i++ ) {
        wait1 = ThreadPool::createWork( check_inc, 1 );
        wait1->add_dependency( id1 );
        tpool.add_work( wait1 );
    }
    tpool.wait_pool_finished();
    ids.clear();
    for ( int i = 0; i < 5; i++ )
        ids.push_back( TPOOL_ADD_WORK( &tpool, sleep_inc2, ( 0.5 ) ) );
    sleep_inc2( 0.002 );
    ThreadPool::WorkItem *work = ThreadPool::createWork( waste_cpu, 100 );
    work->add_dependencies( ids );
    id = tpool.add_work( work, 10 );
    tpool.wait( id );
}


/******************************************************************
 * Test FIFO behavior                                              *
 ******************************************************************/
void test_FIFO( UnitTest &ut, ThreadPool &tpool )
{
    AMP::AMP_MPI comm( MPI_COMM_WORLD );
    int rank    = comm.getRank();
    int size    = comm.getSize();
    const int N = 4000;
    for ( int r = 0; r < size; r++ ) {
        comm.barrier();
        if ( r != rank )
            continue;
        std::vector<ThreadPoolID> ids;
        ids.reserve( N );
        for ( size_t i = 0; i < N; i++ )
            ids.emplace_back( TPOOL_ADD_WORK( &tpool, sleep_inc2, ( 0.001 ) ) );
        bool pass = true;
        while ( tpool.N_queued() > 0 ) {
            int i1 = -1, i2 = ids.size();
            for ( int i = N - 1; i >= 0; i-- ) {
                bool started = ids[i].started();
                if ( started )
                    i1 = std::max<int>( i1, i ); // Last index to processing item
                else
                    i2 = std::min<int>( i2, i ); // First index to queued item
            }
            int diff = i1 == -1 ? 0 : ( i2 - i1 - 1 );
            if ( abs( diff ) > 4 ) {
                printf( "%i %i %i\n", i1, i2, diff );
                pass = pass && abs( i2 - i1 - 1 ) <= 2;
            }
        }
        ids.clear();
        tpool.wait_pool_finished();
        if ( pass )
            ut.passes( "Thread pool behaves as FIFO" );
        else
            ut.failure( "Thread pool does not behave as FIFO" );
    }
}


/******************************************************************
 * Test process affinities                                         *
 ******************************************************************/
void testProcessAffinity( UnitTest &ut )
{
    // Get the number of processors availible
    int N_procs = ThreadPool::getNumberOfProcessors();
    if ( N_procs > 0 )
        ut.passes( "getNumberOfProcessors" );
    else
        ut.failure( "getNumberOfProcessors" );
    printp( "%i processors availible\n", N_procs );

    // Get the processor affinities for the process
    auto cpus = ThreadPool::getProcessAffinity();
    printp( "%i cpus for current process: ", (int) cpus.size() );
    for ( int cpu : cpus )
        printp( "%i ", cpu );
    printp( "\n" );
    if ( !cpus.empty() ) {
        ut.passes( "getProcessAffinity" );
    } else {
#ifdef __APPLE__
        ut.expected_failure( "getProcessAffinity" );
#else
        ut.failure( "getProcessAffinity" );
#endif
    }

    // Test setting the process affinities
    bool pass = false;
    if ( !cpus.empty() && N_procs > 0 ) {
        if ( cpus.size() == 1 ) {
            cpus.resize( N_procs );
            for ( int i = 0; i < N_procs; i++ )
                cpus.push_back( i );
            try {
                ThreadPool::setProcessAffinity( cpus );
            } catch ( ... ) {
            }
            cpus = ThreadPool::getProcessAffinity();
            printp( "%i cpus for current process (updated): ", (int) cpus.size() );
            for ( int cpu : cpus )
                printp( "%i ", cpu );
            printp( "\n" );
            pass = cpus.size() > 1;
        } else {
            auto cpus_orig = cpus;
            try {
                std::vector<int> cpus_tmp( 1, cpus[0] );
                ThreadPool::setProcessAffinity( cpus_tmp );
            } catch ( ... ) {
            }
            cpus = ThreadPool::getProcessAffinity();
            if ( cpus.size() == 1 )
                pass = true;
            try {
                ThreadPool::setProcessAffinity( cpus_orig );
            } catch ( ... ) {
            }
            cpus = ThreadPool::getProcessAffinity();
            if ( cpus.size() != cpus_orig.size() )
                pass = false;
        }
    }
    if ( pass ) {
        ut.passes( "setProcessAffinity" );
    } else {
#ifdef __APPLE__
        ut.expected_failure( "setProcessAffinity" );
#else
        ut.failure( "setProcessAffinity" );
#endif
    }
}


/******************************************************************
 * Test thread affinities                                          *
 ******************************************************************/
void testThreadAffinity( ThreadPool &tpool, UnitTest &ut )
{
    AMP::AMP_MPI comm( MPI_COMM_WORLD );
    int N_threads = tpool.getNumThreads();
    auto cpus     = ThreadPool::getProcessAffinity();
    int rank      = comm.getRank();

    // Test setting the thread affinities
    comm.barrier();
    if ( cpus.size() > 1 ) {
        printp( "Testing thread affinities\n" );
        sleep_ms( 50 );
        // First make sure we can get the thread affinities
        auto procs = ThreadPool::getThreadAffinity();
        if ( procs == cpus ) {
            ut.passes( "getThreadAffinity() matches procs" );
        } else {
            auto msg = stringf( "getThreadAffinity() does not match procs (%i,%i)",
                                static_cast<int>( procs.size() ),
                                static_cast<int>( cpus.size() ) );
            ut.failure( msg );
        }
        auto pass = true;
        for ( int i = 0; i < N_threads; i++ ) {
            auto procs_thread = tpool.getThreadAffinity( i );
            if ( procs_thread != procs ) {
                printp( "%i: Initial thread affinity: ", rank );
                for ( int p : procs_thread )
                    printp( "%i ", p );
                printp( "\n" );
                pass = false;
            }
        }
        if ( pass )
            ut.passes( "getThreadAffinity(thread) matches procs" );
        else
            ut.failure( "getThreadAffinity(thread) does not match procs" );
        // Try to set the thread affinities
        pass = true;
        if ( !procs.empty() ) {
            int N_procs_thread = std::max<int>( (int) cpus.size() / N_threads, 1 );
            for ( int i = 0; i < N_threads; i++ ) {
                std::vector<int> procs_thread( N_procs_thread, -1 );
                for ( int j = 0; j < N_procs_thread; j++ )
                    procs_thread[j] = procs[( i * N_procs_thread + j ) % procs.size()];
                tpool.setThreadAffinity( i, procs_thread );
                sleep_ms( 10 ); // Give time for OS to update thread affinities
                auto procs_thread2 = tpool.getThreadAffinity( i );
                if ( procs_thread2 != procs_thread ) {
                    printp( "%i: Final thread affinity: ", rank );
                    for ( int p : procs_thread )
                        printp( "%i ", p );
                    printp( "\n" );
                    pass = false;
                }
            }
        }
        if ( pass )
            ut.passes( "setThreadAffinity passes" );
        else
            ut.failure( "setThreadAffinity failed to change affinity" );
    }
}


/******************************************************************
 * Test ThreadPool performance                                     *
 ******************************************************************/
void testThreadPoolPerformance( ThreadPool &tpool )
{
    constexpr int N_work    = 2000; // Number of work items
    constexpr int N_it      = 10;   // Number of cycles to run
    constexpr int N_problem = 5;    // Problem size

    const int N_threads = tpool.getNumThreads();
    printp( "\nTesting ThreadPool performance with %i threads:\n", N_threads );

    // Initialize the data
    std::vector<int> data1( N_work, 0 );
    std::vector<int> priority( N_work, 0 );
    std::vector<ThreadPoolID> ids( N_work );
    std::vector<ThreadPool::WorkItem *> work( N_work );
    for ( int i = 0; i < N_work; i++ ) {
        data1[i]    = N_problem;
        priority[i] = i % 128;
    }

    // Test basic cost
    auto start = std::chrono::high_resolution_clock::now();
    for ( int n = 0; n < N_it; n++ ) {
        for ( int i = 0; i < N_work; i++ )
            waste_cpu( data1[i] );
    }
    auto stop     = std::chrono::high_resolution_clock::now();
    int64_t time0 = to_ns( stop - start ) / ( N_it * N_work );
    printp( "   Time for serial item = %i ns\n", time0 );

    // Test the timing creating and running a work item
    printp( "   Testing timmings (creating/running work item):\n" );
    PROFILE_START( "Create/Run work item" );
    tpool.wait_pool_finished();
    int64_t time_create = 0;
    int64_t time_run    = 0;
    int64_t time_delete = 0;
    for ( int n = 0; n < N_it; n++ ) {
        auto t1 = std::chrono::high_resolution_clock::now();
        for ( int i = 0; i < N_work; i++ )
            work[i] = ThreadPool::createWork<void, int>( waste_cpu, data1[i] );
        auto t2 = std::chrono::high_resolution_clock::now();
        for ( int i = 0; i < N_work; i++ )
            work[i]->run();
        auto t3 = std::chrono::high_resolution_clock::now();
        for ( int i = 0; i < N_work; i++ )
            delete work[i];
        auto t4 = std::chrono::high_resolution_clock::now();
        time_create += to_ns( t2 - t1 );
        time_run += to_ns( t3 - t2 );
        time_delete += to_ns( t4 - t3 );
    }
    time_create /= ( N_it * N_work );
    time_run /= ( N_it * N_work );
    time_delete /= ( N_it * N_work );
    PROFILE_STOP( "Create/Run work item" );
    printp( "      create = %i ns\n", time_create );
    printp( "      run    = %i ns\n", time_run );
    printp( "      delete = %i ns\n", time_delete );

    // Test the timing adding a single item
    printp( "   Testing timmings (adding a single item):\n" );
    auto timer_name = Utilities::stringf( "Add single item to tpool (%i threads)", N_threads );
    NULL_USE( timer_name );
    PROFILE_START( timer_name );
    int64_t time_add_single  = 0;
    int64_t time_wait_single = 0;
    for ( int n = 0; n < N_it; n++ ) {
        auto t1 = std::chrono::high_resolution_clock::now();
        for ( int i = 0; i < N_work; i++ )
            ids[i] = TPOOL_ADD_WORK( &tpool, waste_cpu, ( data1[i] ), priority[i] );
        auto t2 = std::chrono::high_resolution_clock::now();
        tpool.wait_all( ids );
        auto t3 = std::chrono::high_resolution_clock::now();
        time_add_single += to_ns( t2 - t1 );
        time_wait_single += to_ns( t3 - t2 );
    }
    time_add_single /= ( N_it * N_work );
    time_wait_single /= ( N_it * N_work );
    PROFILE_STOP( timer_name );
    printp( "      create and add = %i ns\n", time_add_single );
    printp( "      wait = %i ns\n", time_wait_single );

    // Test the timing pre-creating the work items and adding multiple at a time
    printp( "   Testing timmings (adding a block of items):\n" );
    timer_name = Utilities::stringf( "Add block of items to tpool (%i threads)", N_threads );
    NULL_USE( timer_name );
    PROFILE_START( timer_name );
    int64_t time_create_multiple = 0;
    int64_t time_add_multiple    = 0;
    int64_t time_wait_multiple   = 0;
    for ( int n = 0; n < N_it; n++ ) {
        auto t1 = std::chrono::high_resolution_clock::now();
        for ( int i = 0; i < N_work; i++ )
            work[i] = ThreadPool::createWork<void, int>( waste_cpu, data1[i] );
        auto t2   = std::chrono::high_resolution_clock::now();
        auto ids2 = tpool.add_work( work, priority );
        auto t3   = std::chrono::high_resolution_clock::now();
        tpool.wait_all( ids2 );
        auto t4 = std::chrono::high_resolution_clock::now();
        time_create_multiple += to_ns( t2 - t1 );
        time_add_multiple += to_ns( t3 - t2 );
        time_wait_multiple += to_ns( t4 - t3 );
    }
    time_create_multiple /= ( N_it * N_work );
    time_add_multiple /= ( N_it * N_work );
    time_wait_multiple /= ( N_it * N_work );
    PROFILE_STOP( timer_name );
    printp( "      create = %i ns\n", time_create_multiple );
    printp( "      add = %i ns\n", time_add_multiple );
    printp( "      wait = %i ns\n", time_wait_multiple );

    // Estimate the overheads
    int Nt           = std::max( N_threads, 1 );
    int create       = time_create;
    int run          = std::max<int>( time_run - time0, 0 );
    int add_single   = time_add_single - time_create;
    int add_multiple = time_add_multiple;
    int total_single = time_add_single + time_wait_single - ( time0 / Nt );
    int total_multiple =
        time_create_multiple + time_add_multiple + time_wait_multiple - ( time0 / Nt );
    if ( tpool.getNumThreads() == 0 ) {
        add_single -= time0;
        add_multiple -= time0;
    }
    printp( "   Overhead:\n" );
    printp( "      WorkItem create = %i ns\n", create );
    printp( "      WorkItem run    = %i ns\n", run );
    printp( "      Add: serial     = %i ns\n", add_single );
    printp( "      Add: multiple   = %i ns\n", add_multiple );
    printp( "      Total: serial   = %i ns\n", total_single );
    printp( "      Total: multiple = %i ns\n", total_multiple );
}


/******************************************************************
 * The main program                                                *
 ******************************************************************/
void run_tests( UnitTest &ut )
{
    constexpr int N_threads = 4; // Number of threads

    // Disable OS specific warnings for all non-root ranks
    AMP::AMP_MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    if ( rank > 0 )
        ThreadPool::set_OS_warnings( 1 );

    // Test the atomics
    if ( test_atomics() )
        ut.passes( "Atomics passed" );
    else
        ut.failure( "Atomics failed" );

    // Print the size of the thread pool class
    printp( "Size of ThreadPool = %i\n", (int) sizeof( ThreadPool ) );

    // Test process affinities
    comm.barrier();
    testProcessAffinity( ut );
    int N_procs      = ThreadPool::getNumberOfProcessors();
    int N_procs_used = std::min<int>( N_procs, N_threads );
    printp( "%i processors used\n", N_procs_used );

    // Create the thread pool
    comm.barrier();
    printp( "Creating thread pool\n" );
    ThreadPool tpool;
    {
        auto id = TPOOL_ADD_WORK( &tpool, waste_cpu, ( 5 ) );
        if ( id == ThreadPoolID() || !tpool.isValid( id ) )
            ut.failure( "Errors with id" );
    }
    tpool.setNumThreads( N_threads );
    if ( tpool.getNumThreads() == N_threads )
        ut.passes( "Created thread pool" );
    else
        ut.failure( "Failed to create tpool with desired number of threads" );


    // Test setting the thread affinities
    testThreadAffinity( tpool, ut );

    // Print the current processors by thread id
    comm.barrier();
    ThreadPool::set_OS_warnings( 1 );
    print_processor( &tpool );
    launchAndTime( tpool, N_threads, print_processor, &tpool );

    // Test the move operator for thread_id
    test_move_thread_id( ut, tpool );

    // Test calling functions with different number of arguments
    comm.barrier();
    test_function_arguments( &tpool, ut );

    // Check that threads sleep in parallel (does not depend on the number of processors)
    comm.barrier();
    test_sleep_parallel( ut, tpool );

    // Check that the threads are actually working in parallel
    comm.barrier();
    test_work_parallel( ut, tpool );

    // Test first-in-first-out scheduler (also ensures priorities)
    test_FIFO( ut, tpool );

    // Test adding a work item with a dependency
    comm.barrier();
    test_work_dependency( ut, tpool );

    // Run some performance tests
    if ( rank == 0 ) {
        ThreadPool tpool0;
        testThreadPoolPerformance( tpool0 );
        testThreadPoolPerformance( tpool );
    }

    // Run a dependency test that tests a simple case that should keep the thread pool busy
    // Note: Checking the results requires looking at the trace data
    tpool.wait_pool_finished();
    PROFILE_START( "Dependency test" );
    for ( int i = 0; i < 10; i++ ) {
        char msg[3][100];
        snprintf( msg[0], 100, "Item %i-%i", i, 0 );
        snprintf( msg[1], 100, "Item %i-%i", i, 1 );
        snprintf( msg[2], 100, "Item %i-%i", i, 2 );
        auto work  = ThreadPool::createWork( sleep_msg, 0.5, msg[0] );
        auto work1 = ThreadPool::createWork( sleep_msg, 0.1, msg[1] );
        auto work2 = ThreadPool::createWork( sleep_msg, 0.1, msg[2] );
        auto id    = tpool.add_work( work );
        work1->add_dependency( id );
        work2->add_dependency( id );
        tpool.add_work( work1 );
        tpool.add_work( work2 );
    }
    tpool.wait_pool_finished();
    PROFILE_STOP( "Dependency test" );

    // Close the thread pool
    tpool.setNumThreads( 0 );

    // Save the profiling results
    PROFILE_SAVE( "test_thread_pool" );
    PROFILE_DISABLE();

    // Test creating/destroying a thread pool using new
    comm.barrier();
    auto pass = !ThreadPool::is_valid( nullptr );
    try {
        auto tpool2 = new ThreadPool( ThreadPool::MAX_THREADS - 1 );
        if ( tpool2->getNumThreads() != ThreadPool::MAX_THREADS - 1 )
            pass = false;
        if ( !ThreadPool::is_valid( tpool2 ) )
            pass = false;
        delete tpool2;
    } catch ( ... ) {
        pass = false;
    }
    if ( pass )
        ut.passes( "Created/destroyed thread pool with new" );
    else
        ut.failure( "Created/destroyed thread pool with new" );
}
int main( int argc, char *argv[] )
{
    // Initialize MPI and profiler
    AMP::AMPManager::startup( argc, argv );
    UnitTest ut;
    PROFILE_ENABLE( 3 );
    PROFILE_ENABLE_TRACE();
    PROFILE_DISABLE_MEMORY();

    // Run the tests
    run_tests( ut );

    // Shutdown
    ut.report();
    int N_errors = ut.NumFailGlobal();
    ut.reset();
    AMP::AMPManager::shutdown();
    return N_errors;
}
