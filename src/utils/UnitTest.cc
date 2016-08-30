#include "UnitTest.h"
#include "Utilities.h"
#include "AMPManager.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>


namespace AMP {


/********************************************************************
*  Constructor/Destructor                                           *
********************************************************************/
UnitTest::UnitTest() : d_verbose( false )
{
    if ( !AMP::AMPManager::isInitialized() )
        AMP_ERROR( "AMPManager must be initialized first" );
    comm = MPI_COMM_WORLD;
}
UnitTest::~UnitTest() { reset(); }
void UnitTest::reset()
{
    pass_messages.clear();
    fail_messages.clear();
    expected_fail_messages.clear();
}


/********************************************************************
*  Print a global report                                            *
*  Note: only rank 0 will print, all messages will be aggregated    *
********************************************************************/
void UnitTest::report( const int level0 ) const
{
    int size = this->size();
    int rank = this->rank();
    // Give all processors a chance to print any remaining messages
    comm.barrier();
    Utilities::sleepMs( 10 );
    // Broadcast the print level from rank 0
    int level = comm.bcast( level0, 0 );
    if ( level < 0 || level > 2 )
        AMP_ERROR( "Invalid print level" );
    // Perform a global all gather to get the number of failures per processor
    std::vector<int> N_pass( size, 0 );
    std::vector<int> N_fail( size, 0 );
    std::vector<int> N_expected_fail( size, 0 );
    comm.allGather<int>( (int) pass_messages.size(), &N_pass[0] );
    comm.allGather<int>( (int) fail_messages.size(), &N_fail[0] );
    comm.allGather<int>( (int) expected_fail_messages.size(), &N_expected_fail[0] );
    int N_pass_tot          = 0;
    int N_fail_tot          = 0;
    int N_expected_fail_tot = 0;
    for ( int i = 0; i < size; i++ ) {
        N_pass_tot += N_pass[i];
        N_fail_tot += N_fail[i];
        N_expected_fail_tot += N_expected_fail[i];
    }
    NULL_USE( N_fail_tot );
    // Send all messages to rank 0 (if needed)
    std::vector<std::vector<std::string>> pass_messages_rank( comm.getSize() );
    std::vector<std::vector<std::string>> fail_messages_rank( comm.getSize() );
    std::vector<std::vector<std::string>> expected_fail_rank( comm.getSize() );
    // Get the pass messages
    if ( ( level == 1 && N_pass_tot <= 20 ) || level == 2 ) {
        if ( rank == 0 ) {
            // Rank 0 should receive all messages
            for ( int i = 0; i < size; i++ ) {
                if ( i == 0 )
                    pass_messages_rank[i] = pass_messages;
                else if ( N_pass[i] > 0 )
                    pass_messages_rank[i] = unpack_message_stream( i, 1 );
            }
        } else if ( !pass_messages.empty() ) {
            // All other ranks send their message (use non-blocking communication)
            pack_message_stream( pass_messages, 0, 1 );
        }
    }
    // Get the fail messages
    if ( level == 1 || level == 2 ) {
        if ( rank == 0 ) {
            // Rank 0 should receive all messages
            for ( int i = 0; i < size; i++ ) {
                if ( i == 0 )
                    fail_messages_rank[i] = fail_messages;
                else if ( N_fail[i] > 0 )
                    fail_messages_rank[i] = unpack_message_stream( i, 2 );
            }
        } else if ( !fail_messages.empty() ) {
            // All other ranks send their message (use non-blocking communication)
            pack_message_stream( fail_messages, 0, 2 );
        }
    }
    // Get the expected_fail messages
    if ( ( level == 1 && N_expected_fail_tot <= 50 ) || level == 2 ) {
        if ( rank == 0 ) {
            // Rank 0 should receive all messages
            for ( int i = 0; i < size; i++ ) {
                if ( i == 0 )
                    expected_fail_rank[i] = expected_fail_messages;
                else if ( N_expected_fail[i] > 0 )
                    expected_fail_rank[i] = unpack_message_stream( i, 3 );
            }
        } else if ( !expected_fail_messages.empty() ) {
            // All other ranks send their message (use non-blocking communication)
            pack_message_stream( expected_fail_messages, 0, 3 );
        }
    }
    // Print the results of all messages (only rank 0 will print)
    if ( rank == 0 ) {
        std::cout << std::endl;
        // Print the passed tests
        std::cout << "Tests passed" << std::endl;
        if ( level == 0 || ( level == 1 && N_pass_tot > 20 ) ) {
            // We want to print a summary
            if ( size > 8 ) {
                // Print 1 summary for all processors
                std::cout << "     " << N_pass_tot
                          << " tests passed (use report level 2 for more detail)" << std::endl;
            } else {
                // Print a summary for each processor
                for ( int i = 0; i < size; i++ )
                    std::cout << "     " << N_pass[i] << " tests passed (proc " << i
                              << ") (use report level 2 for more detail)" << std::endl;
            }
        } else {
            // We want to print all messages
            for ( int i = 0; i < size; i++ ) {
                AMP_ASSERT( (int) pass_messages_rank[i].size() == N_pass[i] );
                if ( N_pass[i] > 0 ) {
                    std::cout << "     Proccessor " << i << ":" << std::endl;
                    for ( auto &elem : pass_messages_rank[i] )
                        std::cout << "        " << elem << std::endl;
                }
            }
        }
        std::cout << std::endl;
        // Print the tests that failed
        std::cout << "Tests failed" << std::endl;
        if ( level == 0 ) {
            // We want to print a summary
            if ( size > 8 ) {
                // Print 1 summary for all processors
                std::cout << "     " << N_pass_tot
                          << " tests failed (use report level 2 for more detail)" << std::endl;
            } else {
                // Print a summary for each processor
                for ( int i = 0; i < size; i++ )
                    std::cout << "     " << N_fail[i] << " tests failed (proc " << i
                              << ") (use report level 1 or 2 for more detail)" << std::endl;
            }
        } else {
            // We want to print all messages
            for ( int i = 0; i < size; i++ ) {
                AMP_ASSERT( (int) fail_messages_rank[i].size() == N_fail[i] );
                if ( N_fail[i] > 0 ) {
                    std::cout << "     Processor " << i << ":" << std::endl;
                    for ( auto &elem : fail_messages_rank[i] )
                        std::cout << "        " << elem << std::endl;
                }
            }
        }
        std::cout << std::endl;
        // Print the tests that expected failed
        std::cout << "Tests expected failed" << std::endl;
        if ( level == 0 || ( level == 1 && N_expected_fail_tot > 50 ) ) {
            // We want to print a summary
            if ( size > 8 ) {
                // Print 1 summary for all processors
                std::cout << "     " << N_expected_fail_tot
                          << " tests expected failed (use report level 2 for more detail)"
                          << std::endl;
            } else {
                // Print a summary for each processor
                for ( int i = 0; i < size; i++ )
                    std::cout << "     " << N_expected_fail[i] << " tests expected failed (proc "
                              << i << ") (use report level 1 or 2 for more detail)" << std::endl;
            }
        } else {
            // We want to print all messages
            for ( int i = 0; i < size; i++ ) {
                AMP_ASSERT( (int) expected_fail_rank[i].size() == N_expected_fail[i] );
                if ( N_expected_fail[i] > 0 ) {
                    std::cout << "     Processor " << i << ":" << std::endl;
                    for ( auto &elem : expected_fail_rank[i] )
                        std::cout << "        " << elem << std::endl;
                }
            }
        }
        std::cout << std::endl;
    }
    // Add a barrier to synchronize all processors (rank 0 is much slower)
    comm.barrier();
    Utilities::sleepMs( 10 ); // Need a brief pause to allow any printing to finish
}


/********************************************************************
*  Pack and send the given messages                                 *
********************************************************************/
void UnitTest::pack_message_stream( const std::vector<std::string> &messages,
                                    const int rank,
                                    const int tag ) const
{
    // Get the size of the messages
    int N_messages   = (int) messages.size();
    int *msg_size    = new int[N_messages];
    int msg_size_tot = 0;
    for ( int i = 0; i < N_messages; i++ ) {
        msg_size[i] = (int) messages[i].size();
        msg_size_tot += msg_size[i];
    }
    // Allocate space for the message stream
    int size_data = ( N_messages + 1 ) * sizeof( int ) + msg_size_tot;
    char *data    = new char[size_data];
    // Pack the message stream
    int *tmp = (int *) data;
    tmp[0]   = N_messages;
    for ( int i    = 0; i < N_messages; i++ )
        tmp[i + 1] = msg_size[i];
    int k          = ( N_messages + 1 ) * sizeof( int );
    for ( int i = 0; i < N_messages; i++ ) {
        messages[i].copy( &data[k], msg_size[i] );
        k += msg_size[i];
    }
    // Send the message stream (using a non-blocking send)
    MPI_Request request = comm.Isend( data, size_data, rank, tag );
    // Wait for the communication to send and free the temporary memory
    AMP::AMP_MPI::wait( request );
    delete[] data;
    delete[] msg_size;
}


/********************************************************************
*  Receive and unpack a message stream                              *
********************************************************************/
std::vector<std::string> UnitTest::unpack_message_stream( const int rank, const int tag ) const
{
    // Probe the message to get the message size
    int size_data = comm.probe( rank, tag );
    // Allocate memory to receive the data
    char *data = new char[size_data];
    // Receive the data (using a non-blocking receive)
    MPI_Request request = comm.Irecv( data, size_data, rank, tag );
    // Wait for the communication to be received
    AMP::AMP_MPI::wait( request );
    // Unpack the message stream
    int *tmp       = (int *) data;
    int N_messages = tmp[0];
    int *msg_size  = &tmp[1];
    std::vector<std::string> messages( N_messages );
    int k = ( N_messages + 1 ) * sizeof( int );
    for ( int i = 0; i < N_messages; i++ ) {
        messages[i] = std::string( &data[k], msg_size[i] );
        k += msg_size[i];
    }
    // Delete the temporary memory
    delete[] data;
    return messages;
}


/********************************************************************
*  Other functions                                                  *
********************************************************************/
int UnitTest::rank() const { return comm.getRank(); }
int UnitTest::size() const { return comm.getSize(); }
size_t UnitTest::NumPassGlobal() const { return comm.sumReduce<size_t>( NumPassLocal() ); }
size_t UnitTest::NumFailGlobal() const { return comm.sumReduce<size_t>( NumFailLocal() ); }
size_t UnitTest::NumExpectedFailGlobal() const
{
    return comm.sumReduce<size_t>( NumExpectedFailLocal() );
}


} // AMP namespace
