#include "utils/LapackWrappers.h"
#include "utils/Utilities.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace AMP;


// The main function
int main( int, char *[] )
{
    int N_errors = 0;

    // Print the machine specifics
    printf( "\nDouble precision machine parameters\n" );
    Lapack<double>::print_machine_parameters();
    printf( "\nSingle precision machine parameters\n" );
    Lapack<float>::print_machine_parameters();

    // Run the basic tests
    printf( "\nRunning double precision basic tests\n" );
    N_errors += Lapack<double>::run_all_test();
    if ( N_errors == 0 ) {
        printf( "  passed dp tests\n" );
    } else {
        printf( "failed %d dp tests\n", N_errors );
    }
    printf( "\nRunning single precision basic tests\n" );
    N_errors += Lapack<float>::run_all_test();
    if ( N_errors == 0 ) {
        printf( "  passed sp tests\n" );
    } else {
        printf( "failed %d sp tests \n", N_errors );
    }

    // Get the times for the tests
    printf( "\nGetting double precision test times\n" );
    const char *dptests[] = { "dcopy",  "dscal",  "dnrm2",  "dasum",  "ddot",   "daxpy",
                              "dgemv",  "dgemm",  "dgesv",  "dgtsv",  "dgbsv",  "dgetrf",
                              "dgttrf", "dgbtrf", "dgetrs", "dgttrs", "dgbtrs", "dgetri" };
    const int dpN[] = { 500, 500, 500, 500, 500, 100, 100, 100, 100,
                        500, 500, 100, 500, 500, 100, 500, 500, 100 };
    for ( size_t i = 0; i < sizeof( dptests ) / sizeof( char * ); i++ ) {
        double t1    = Utilities::time();
        double error = 0;
        int err      = Lapack<double>::run_test( dptests[i], dpN[i], error );
        double t2    = Utilities::time();
        int us       = static_cast<int>( 1e6 * ( t2 - t1 ) / dpN[i] );
        printf(
            "%7s:  %s:  %5i us  (%e)\n", dptests[i], err == 0 ? "passed" : "failed", us, error );
        N_errors += err;
    }
    if ( N_errors == 0 ) {
        printf( "  passed dp timing tests\n" );
    } else {
        printf( "failed %d dp timing tests\n", N_errors );
    }

    printf( "\nGetting single precision test times\n" );
    const char *sptests[] = { "scopy",  "sscal",  "snrm2",  "sasum",  "sdot",   "saxpy",
                              "sgemv",  "sgemm",  "sgesv",  "sgtsv",  "sgbsv",  "sgetrf",
                              "sgttrf", "sgbtrf", "sgetrs", "sgttrs", "sgbtrs", "sgetri" };
    const int spN[] = { 500, 500, 500, 500, 500, 100, 100, 100, 100,
                        500, 500, 100, 500, 500, 100, 500, 500, 100 };
    for ( size_t i = 0; i < sizeof( sptests ) / sizeof( char * ); i++ ) {
        double t1   = Utilities::time();
        float error = 0;
        int err     = Lapack<float>::run_test( sptests[i], spN[i], error );
        double t2   = Utilities::time();
        int us      = static_cast<int>( 1e6 * ( t2 - t1 ) / spN[i] );
        printf(
            "%7s:  %s:  %5i us  (%e)\n", sptests[i], err == 0 ? "passed" : "failed", us, error );
        N_errors += err;
    }
    if ( N_errors == 0 ) {
        printf( "  passed sp timing tests\n" );
    } else {
        printf( "failed %d sp timing tests\n", N_errors );
    }

    // Finished
    if ( N_errors == 0 )
        std::cout << "\nAll tests passed\n";
    else
        std::cout << "\nSome tests failed\n";
    return N_errors == 0 ? 0 : 1;
}
