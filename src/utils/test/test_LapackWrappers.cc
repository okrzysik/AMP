#include "LapackWrappers.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
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


// The main function
int main( int, char *[] )
{
    // Number of times to run the tests for timing results
    // Note: the tests are designed to take ~ the same time/test
    const int N_test = 50;

    // Store the number of errors
    int N_errors = 0;

    // Set the number of threads
    Lapack<double>::set_num_threads( 1 );

    // Print the lapack version
    std::cout << Lapack<double>::info();

    // Run the basic tests
    printf( "\nRunning double precision basic tests\n" );
    int err0 = Lapack<double>::run_all_test();
    if ( err0 == 0 ) {
        printf( "  passed dp tests\n" );
    } else {
        printf( "  failed %d dp tests\n", err0 );
        N_errors += err0;
    }
    printf( "\nRunning single precision basic tests\n" );
    err0 = Lapack<float>::run_all_test();
    if ( err0 == 0 ) {
        printf( "  passed sp tests\n" );
    } else {
        printf( "  failed %d sp tests \n", err0 );
        N_errors += err0;
    }

    // Get the times for the tests (double)
    {
        printf( "\nGetting double precision test times\n" );
        auto tests = Lapack<double>::list_all_tests();
        std::vector<double> error( tests.size() );
        std::vector<int> time( tests.size() );
        std::vector<int> err( tests.size() );
        int N_err = 0;
        for ( size_t i = 0; i < tests.size(); i++ ) {
            run_test<double>( tests[i], N_test, time[i], error[i], err[i] );
            printf( "%7s:  %s:  %5i us  (%e)\n",
                    tests[i].c_str(),
                    err[i] == 0 ? "passed" : "failed",
                    time[i],
                    error[i] );
            N_err += err[i];
        }
        if ( N_err == 0 ) {
            printf( "  passed dp timing tests\n" );
        } else {
            printf( "  failed %d dp timing tests\n", N_err );
            N_errors += N_err;
        }
    }

    // Get the times for the tests (single)
    {
        printf( "\nGetting single precision test times\n" );
        auto tests = Lapack<float>::list_all_tests();
        std::vector<float> error( tests.size() );
        std::vector<int> time( tests.size() );
        std::vector<int> err( tests.size() );
        int N_err = 0;
        for ( size_t i = 0; i < tests.size(); i++ ) {
            run_test<float>( tests[i], N_test, time[i], error[i], err[i] );
            printf( "%7s:  %s:  %5i us  (%e)\n",
                    tests[i].c_str(),
                    err[i] == 0 ? "passed" : "failed",
                    time[i],
                    error[i] );
            N_err += err[i];
        }
        if ( N_err == 0 ) {
            printf( "  passed sp timing tests\n" );
        } else {
            printf( "  failed %d sp timing tests\n", N_err );
            N_errors += N_err;
        }
    }

    // Run the tests in parallel to check for parallel bugs
    {
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
                threads[j] = std::thread( run_test<double>,
                                          test,
                                          N_test,
                                          std::ref( time2[j] ),
                                          std::ref( error2[j] ),
                                          std::ref( err2[j] ) );
            for ( int j = 0; j < N_threads; j++ )
                threads[j].join();
            auto t2   = std::chrono::system_clock::now();
            bool pass = true;
            for ( int j = 0; j < N_threads; j++ )
                pass = pass && err2[j] == 0;
            auto us = static_cast<int>( diff( t1, t2 ) / ( N_test * N_threads ) );
            printf( "%7s:  %s:  %5i us\n", test.c_str(), pass ? "passed" : "failed", us );
            N_errors += ( pass ? 0 : 1 );
        }
    }

    // Finished
    if ( N_errors == 0 )
        std::cout << "\nAll tests passed\n";
    else
        std::cout << "\nSome tests failed\n";
    return N_errors == 0 ? 0 : 1;
}
