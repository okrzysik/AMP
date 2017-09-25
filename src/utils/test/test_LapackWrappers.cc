#include "LapackWrappers.h"

#include <chrono>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <thread>


// Call Lapack::run_test
template<typename TYPE>
void run_test( const char *routine, int N, int &N_errors, double &error )
{
    N_errors = Lapack<TYPE>::run_test( routine, N, error );
}


// Get the time difference in us
static inline int64_t diff( std::chrono::time_point<std::chrono::system_clock> t1,
                            std::chrono::time_point<std::chrono::system_clock> t2 )
{
    return static_cast<int64_t>( 1e6 * std::chrono::duration<double>( t2 - t1 ).count() );
}


// The main function
int main( int, char *[] )
{
    int N_errors = 0;

    // Print the lapack version
    Lapack<double>::print_lapack_version();
    printf( "\n" );

    // Print the machine specifics
    printf( "\nDouble precision machine parameters\n" );
    Lapack<double>::print_machine_parameters();
    printf( "\nSingle precision machine parameters\n" );
    Lapack<float>::print_machine_parameters();

    // Run the basic tests
    printf( "\nRunning double precision basic tests\n" );
    int err = Lapack<double>::run_all_test();
    if ( err == 0 ) {
        printf( "  passed dp tests\n" );
    } else {
        printf( "  failed %d dp tests\n", err );
        N_errors += err;
    }
    printf( "\nRunning single precision basic tests\n" );
    err = Lapack<float>::run_all_test();
    if ( err == 0 ) {
        printf( "  passed sp tests\n" );
    } else {
        printf( "  failed %d sp tests \n", err );
        N_errors += err;
    }

    // Get the times for the tests (double)
    printf( "\nGetting double precision test times\n" );
    const char *dptests[] = { "dcopy",  "dscal",  "dnrm2",  "dasum",  "ddot",   "daxpy",
                              "dgemv",  "dgemm",  "dgesv",  "dgtsv",  "dgbsv",  "dgetrf",
                              "dgttrf", "dgbtrf", "dgetrs", "dgttrs", "dgbtrs", "dgetri" };
    const int dpN[]       = { 500, 500, 500, 500, 500, 100, 100, 100, 100,
                        500, 500, 100, 500, 500, 100, 500, 500, 100 };
    int N_err             = 0;
    for ( size_t i = 0; i < sizeof( dptests ) / sizeof( char * ); i++ ) {
        auto t1 = std::chrono::system_clock::now();
        double error;
        int err = Lapack<double>::run_test( dptests[i], dpN[i], error );
        auto t2 = std::chrono::system_clock::now();
        int us  = static_cast<int>( diff( t1, t2 ) / dpN[i] );
        printf(
            "%7s:  %s:  %5i us  (%e)\n", dptests[i], err == 0 ? "passed" : "failed", us, error );
        N_err += err;
    }
    if ( N_err == 0 ) {
        printf( "  passed dp timing tests\n" );
    } else {
        printf( "  failed %d dp timing tests\n", N_err );
        N_errors += N_err;
    }

    // Get the times for the tests (single)
    printf( "\nGetting single precision test times\n" );
    const char *sptests[] = { "scopy",  "sscal",  "snrm2",  "sasum",  "sdot",   "saxpy",
                              "sgemv",  "sgemm",  "sgesv",  "sgtsv",  "sgbsv",  "sgetrf",
                              "sgttrf", "sgbtrf", "sgetrs", "sgttrs", "sgbtrs", "sgetri" };
    const int spN[]       = { 500, 500, 500, 500, 500, 100, 100, 100, 100,
                        500, 500, 100, 500, 500, 100, 500, 500, 100 };
    N_err                 = 0;
    for ( size_t i = 0; i < sizeof( sptests ) / sizeof( char * ); i++ ) {
        auto t1 = std::chrono::system_clock::now();
        float error;
        int err = Lapack<float>::run_test( sptests[i], spN[i], error );
        auto t2 = std::chrono::system_clock::now();
        int us  = static_cast<int>( diff( t1, t2 ) / spN[i] );
        printf(
            "%7s:  %s:  %5i us  (%e)\n", sptests[i], err == 0 ? "passed" : "failed", us, error );
        N_err += err;
    }
    if ( N_err == 0 ) {
        printf( "  passed sp timing tests\n" );
    } else {
        printf( "  failed %d sp timing tests\n", N_err );
        N_errors += N_err;
    }

    // Run the tests in parallel to check for parallel bugs
    printf( "\nRunning parallel tests\n" );
    int N_threads = 8;
    for ( size_t i = 0; i < sizeof( dptests ) / sizeof( char * ); i++ ) {
        auto t1 = std::chrono::system_clock::now();
        std::thread threads[128];
        int N_errors_thread[128];
        double error_thread[128];
        for ( int j = 0; j < N_threads; j++ )
            threads[j] = std::thread( run_test<double>,
                                      dptests[i],
                                      dpN[i],
                                      std::ref( N_errors_thread[j] ),
                                      std::ref( error_thread[j] ) );
        for ( int j = 0; j < N_threads; j++ )
            threads[j].join();
        auto t2   = std::chrono::system_clock::now();
        bool pass = true;
        for ( int j = 0; j < N_threads; j++ )
            pass = pass && N_errors_thread[j] == 0;
        int us = static_cast<int>( diff( t1, t2 ) / ( dpN[i] * N_threads ) );
        printf( "%7s:  %s:  %5i us\n", dptests[i], pass ? "passed" : "failed", us );
        N_errors += ( pass ? 0 : 1 );
    }


    // Finished
    if ( N_errors == 0 )
        std::cout << "\nAll tests passed\n";
    else
        std::cout << "\nSome tests failed\n";
    return N_errors == 0 ? 0 : 1;
}
