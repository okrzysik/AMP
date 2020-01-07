// This is a helper test that will mimic the type of thread placement that nrdf uses
// and print the actual affinities of the threads.  This is useful for determining /
// verifying the optimum placement on large machines (titan).

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/threadpool/ThreadPool.h"


using namespace AMP;


/************************************************************************
 * Get/Create the thread pool database and the local thread pool to use  *
 ************************************************************************/
std::shared_ptr<ThreadPool> create_thread_pool( std::shared_ptr<AMP::Database> tpool_db )
{
    std::shared_ptr<ThreadPool> tpool;
    // Check the ThreadPool database and create the ThreadPool
    if ( tpool_db != nullptr ) {
        AMP_ASSERT( tpool_db->keyExists( "N_threads" ) );
        AMP_ASSERT( tpool_db->keyExists( "share_tpool" ) );
        if ( !tpool_db->getScalar<bool>( "share_tpool" ) ) {
            AMP_ASSERT( tpool_db->keyExists( "N_threads_E" ) );
            AMP_ASSERT( tpool_db->keyExists( "N_threads_T" ) );
        }
        // Get the number threads needed and the processors availible
        int N_threads     = tpool_db->getWithDefault( "N_threads", 0 );
        int N_threads_max = std::max( N_threads, 1 );
        if ( !tpool_db->getScalar<bool>( "share_tpool" ) ) {
            int N_threads_T = tpool_db->getWithDefault( "N_threads_T", 0 );
            int N_threads_E = tpool_db->getWithDefault( "N_threads_E", 0 );
            N_threads_max   = std::max( N_threads_max, N_threads_T + N_threads_E );
        }
        // First we need to check and modify the affinities of the current process
        int method = tpool_db->getWithDefault( "Load_balance_method", 1 );
        int N_min  = tpool_db->getWithDefault( "N_min", -1 );
        int N_max  = tpool_db->getWithDefault( "N_max", -1 );
        if ( N_min == -1 )
            N_min = std::min( N_threads_max, N_max );
        auto procs = ThreadPool::getProcessAffinity();
        if ( procs.empty() ) {
            // The current OS does not support idenifying or setting the affinity
            AMP::pout << "Warning: OS does not support process affinities\n";
        } else if ( tpool_db->keyExists( "procs" ) ) {
            int offset = procs[0];
            procs      = tpool_db->getVector<int>( "procs" );
            for ( int &proc : procs )
                proc += offset;
            ThreadPool::setProcessAffinity( procs );
        } else {
            AMP::AMP_MPI::balanceProcesses(
                AMP::AMP_MPI( AMP_COMM_WORLD ), method, std::vector<int>(), N_min, N_max );
        }
        procs = ThreadPool::getProcessAffinity();
        // Create the thread pool
        if ( N_threads == 0 )
            return tpool;
        if ( method == 1 )
            tpool.reset( new ThreadPool( N_threads ) );
        else
            tpool.reset( new ThreadPool( N_threads, "independent", procs ) );
    }
    return tpool;
}


/************************************************************************
 * Create the nested thread pools                                        *
 ************************************************************************/
void create_nested_thread_pool( std::shared_ptr<AMP::Database> tpool_db,
                                std::shared_ptr<ThreadPool> &tpool_E,
                                std::shared_ptr<ThreadPool> &tpool_T )
{
    int N_threads_E = tpool_db->getWithDefault( "N_threads_E", 0 );
    int N_threads_T = tpool_db->getWithDefault( "N_threads_T", 0 );
    if ( tpool_db->getScalar<bool>( "share_tpool" ) || N_threads_E + N_threads_T == 0 ) {
        tpool_E.reset();
        tpool_T.reset();
    } else {
        auto procs = ThreadPool::getProcessAffinity();
        std::vector<int> procs_E, procs_T;
        if ( N_threads_E + N_threads_T == (int) procs.size() ) {
            for ( int i = 0; i < N_threads_E; ++i )
                procs_E.push_back( procs[i] );
            for ( int i = 0; i < N_threads_T; ++i )
                procs_T.push_back( procs[N_threads_E + i] );
        } else if ( N_threads_E + N_threads_T < (int) procs.size() ) {
            int n_E = ( N_threads_E * procs.size() ) / ( N_threads_E + N_threads_T );
            int n_T = procs.size() - n_E;
            for ( int i = 0; i < n_E; ++i )
                procs_E.push_back( procs[i] );
            for ( int i = 0; i < n_T; ++i )
                procs_T.push_back( procs[n_E + i] );
        } else if ( N_threads_E + N_threads_T > (int) procs.size() ) {
            procs_E = procs;
            procs_T = procs;
        } else {
            AMP_ERROR( "Internal error" );
        }
        if ( !procs_E.empty() )
            tpool_E.reset( new ThreadPool( N_threads_E, "independent", procs_E ) );
        if ( !procs_T.empty() )
            tpool_T.reset( new ThreadPool( N_threads_T, "independent", procs_T ) );
    }
}


/******************************************************************
 * Return the thread affinities for a thread pool                  *
 ******************************************************************/
std::vector<std::vector<bool>> get_thread_affinity( std::shared_ptr<ThreadPool> tpool )
{
    std::vector<std::vector<bool>> affinity;
    if ( tpool != nullptr ) {
        int N_procs = tpool->getNumberOfProcessors();
        affinity.resize( tpool->getNumThreads(), std::vector<bool>( N_procs, false ) );
        for ( int i = 0; i < tpool->getNumThreads(); i++ ) {
            auto procs = tpool->getThreadAffinity( i );
            for ( int proc : procs )
                affinity[i][proc] = true;
        }
    }
    return affinity;
}


/******************************************************************
 * Print the affinities                                            *
 ******************************************************************/
inline void print_affinity( const std::vector<bool> &affinity )
{
    for ( bool j : affinity )
        printf( "%4s  ", j ? "X" : "" );
    printf( "\n" );
}


/******************************************************************
 * The main program                                                *
 ******************************************************************/
int main( int argc, char *argv[] )
{

    // Initialize MPI and set the error handlers
    AMP::AMP_MPI::start_MPI( argc, argv );
    AMP::UnitTest ut;

    { // Limit scope

        // Create the tpool database
        auto tpool_db = std::make_shared<AMP::Database>( "ThreadPool" );
        if ( argc == 1 ) {
            // Create a default database
            tpool_db.reset( new AMP::Database( "ThreadPool" ) );
            tpool_db->putScalar( "Load_balance_method", 2 );
            tpool_db->putScalar( "N_min", -1 );
            tpool_db->putScalar( "N_max", -1 );
            tpool_db->putScalar( "N_threads", 2 );
            tpool_db->putScalar( "share_tpool", false );
            tpool_db->putScalar( "N_threads_E", 0 );
            tpool_db->putScalar( "N_threads_T", 0 );
        } else if ( argc == 2 ) {
            // Read the database from the file
            tpool_db = AMP::Database::parseInputFile( argv[1] );
        } else {
            AMP_ERROR( "Invalid number of inputs" );
        }

        // Create the threadpools
        auto tpool = create_thread_pool( tpool_db );
        std::shared_ptr<ThreadPool> tpool_E, tpool_T;
        create_nested_thread_pool( tpool_db, tpool_E, tpool_T );

        // Get the list of affinities for the current process and all threads
        int N_procs = tpool->getNumberOfProcessors();
        auto procs2 = tpool->getProcessAffinity();
        std::vector<bool> affinity_process( N_procs, false );
        for ( int j : procs2 )
            affinity_process[j] = true;
        auto affinity_thread = get_thread_affinity( tpool );
        auto affinity_E      = get_thread_affinity( tpool_E );
        auto affinity_T      = get_thread_affinity( tpool_T );

        // Get subcomms for the different nodes
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        auto nodeComm = globalComm.splitByNode();

        // Get the rank of node root and the set of node roots
        int nodeRoot = nodeComm.minReduce( globalComm.getRank() );
        std::set<int> nodeSet;
        nodeSet.insert( nodeRoot );
        globalComm.setGather( nodeSet );
        int nodeIndex = 0;
        for ( auto it = nodeSet.begin(); it != nodeSet.end(); ++it, ++nodeIndex ) {
            globalComm.barrier(); // Syncronize all processors for output
            if ( *it == nodeRoot ) {
                for ( int i = 0; i < nodeComm.getSize(); i++ ) {
                    globalComm.barrier(); // Syncronize all processors for output
                    if ( i == nodeComm.getRank() ) {
                        // Print the header
                        if ( i == 0 ) {
                            printf( "Node %i: %s\n", nodeIndex, globalComm.getNodeName().c_str() );
                            printf( "%16s", "" );
                            for ( int j = 0; j < N_procs; j++ ) {
                                printf( "%4i  ", j );
                            }
                            printf( "\n" );
                        }
                        // Print the process and thread affinities
                        printf( "Rank %5i, %2i: ", globalComm.getRank(), nodeComm.getRank() );
                        print_affinity( affinity_process );
                        for ( int k = 0; k < (int) affinity_thread.size(); k++ ) {
                            printf( "   thread %2i:   ", k );
                            print_affinity( affinity_thread[k] );
                        }
                        for ( int k = 0; k < (int) affinity_E.size(); k++ ) {
                            printf( "   thread_E %i:  ", k );
                            print_affinity( affinity_E[k] );
                        }
                        for ( int k = 0; k < (int) affinity_T.size(); k++ ) {
                            printf( "   thread_T %i:  ", k );
                            print_affinity( affinity_T[k] );
                        }
                        if ( nodeComm.getRank() == nodeComm.getSize() - 1 )
                            printf( "\n" );
                    }
                }
            }
        }

    } // limit scope

    // Finished
    AMP::AMP_MPI::stop_MPI();
    return 0;
}
