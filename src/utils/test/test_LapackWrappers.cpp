#include "LapackWrappers.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>


// Get the time difference in us
static inline int64_t diff( std::chrono::time_point<std::chrono::system_clock> t1,
    std::chrono::time_point<std::chrono::system_clock> t2 )
{
    return static_cast<int64_t>( 1e6 * std::chrono::duration<double>( t2 - t1 ).count() );
}


// Call Lapack::run_test
template<typename TYPE>
void run_test( const std::string &routine, int N, int &us, TYPE &error, int &err )
{
    auto t1 = std::chrono::system_clock::now();
    err     = Lapack<TYPE>::run_test( routine, N, error );
    auto t2 = std::chrono::system_clock::now();
    us      = static_cast<int>( diff( t1, t2 ) / N );
}


// Print test result
void printTest(
    bool pass, bool print, const std::string &name, int us, double error, int &N_errors )
{
    if ( pass && print ) {
        printf( "%7s:  passed:  %5i us  (%e)\n", name.c_str(), us, error );
    } else if ( !pass ) {
        printf( "%7s:  failed:  %5i us  (%e)\n", name.c_str(), us, error );
        N_errors++;
    }
}


// The main function to run all tests
int runAll( bool print_all )
{
    // Number of times to run the tests for timing results
    // Note: the tests are designed to take ~ the same time/test
    const int N_test = 50;

    // Store the number of errors
    int N_errors = 0;

    // Run the basic tests
    if ( print_all )
        printf( "\nRunning double precision basic tests\n" );
    int err0 = Lapack<double>::run_all_test();
    if ( err0 == 0 && print_all ) {
        printf( "  passed dp tests\n" );
    } else if ( err0 != 0 ) {
        printf( "  failed %d dp tests\n", err0 );
        N_errors += err0;
    }
    if ( print_all )
        printf( "\nRunning single precision basic tests\n" );
    err0 = Lapack<float>::run_all_test();
    if ( err0 == 0 && print_all ) {
        printf( "  passed sp tests\n" );
    } else if ( err0 != 0 ) {
        printf( "  failed %d sp tests \n", err0 );
        N_errors += err0;
    }

    // Get the times for the tests (double)
    {
        if ( print_all )
            printf( "\nGetting double precision test times\n" );
        auto tests = Lapack<double>::list_all_tests();
        std::vector<double> error( tests.size() );
        std::vector<int> time( tests.size() );
        std::vector<int> err( tests.size() );
        int N_err = 0;
        for ( size_t i = 0; i < tests.size(); i++ ) {
            run_test<double>( tests[i], N_test, time[i], error[i], err[i] );
            printTest( err[i] == 0, print_all, tests[i], time[i], error[i], N_err );
        }
        if ( N_err == 0 && print_all ) {
            printf( "  passed dp timing tests\n" );
        } else if ( N_err != 0 ) {
            printf( "  failed %d dp timing tests\n", N_err );
            N_errors += N_err;
        }
    }

    // Get the times for the tests (single)
    {
        if ( print_all )
            printf( "\nGetting single precision test times\n" );
        auto tests = Lapack<float>::list_all_tests();
        std::vector<float> error( tests.size() );
        std::vector<int> time( tests.size() );
        std::vector<int> err( tests.size() );
        int N_err = 0;
        for ( size_t i = 0; i < tests.size(); i++ ) {
            run_test<float>( tests[i], N_test, time[i], error[i], err[i] );
            printTest( err[i] == 0, print_all, tests[i], time[i], error[i], N_err );
        }
        if ( N_err == 0 && print_all ) {
            printf( "  passed sp timing tests\n" );
        } else if ( N_err != 0 ) {
            printf( "  failed %d sp timing tests\n", N_err );
            N_errors += N_err;
        }
    }

    // Run the tests in parallel to check for parallel bugs
    {
        if ( print_all )
            printf( "\nRunning parallel tests\n" );
        int N_threads = 8;
        auto tests    = Lapack<double>::list_all_tests();
        for ( auto &test : tests ) {
            auto t1 = std::chrono::system_clock::now();
            std::thread threads[128];
            int time2[128];
            int err2[128];
            double error2[128];
            for ( int j = 0; j < N_threads; j++ )
                threads[j] = std::thread( run_test<double>, test, N_test, std::ref( time2[j] ),
                    std::ref( error2[j] ), std::ref( err2[j] ) );
            for ( int j = 0; j < N_threads; j++ )
                threads[j].join();
            auto t2      = std::chrono::system_clock::now();
            bool pass    = true;
            double error = 0;
            for ( int j = 0; j < N_threads; j++ ) {
                pass  = pass && err2[j] == 0;
                error = std::max( error, error2[j] );
            }
            auto us = static_cast<int>( diff( t1, t2 ) / ( N_test * N_threads ) );
            printTest( pass, print_all, test, us, error, N_errors );
        }
    }

    return N_errors;
}


// main
int main( int argc, char *argv[] )
{
    // Read inputs
    int N      = 1;
    bool print = true;
    if ( argc == 2 ) {
        print = false;
        N     = std::stoi( argv[1] );
    } else if ( argc > 2 ) {
        std::cerr << "Invald call";
        return 1;
    }

    // Set the number of threads
    Lapack<double>::set_num_threads( 1 );

    // Print the lapack version
    if ( print )
        std::cout << Lapack<double>::info();

    // Run the tests
    int N_errors = 0;
    for ( int i = 0; i < N; i++ )
        N_errors += runAll( print );

    // Finished
    if ( N_errors == 0 )
        std::cout << "\nAll tests passed\n";
    else
        std::cout << "\nSome tests failed\n";
    return N_errors == 0 ? 0 : 1;
}
