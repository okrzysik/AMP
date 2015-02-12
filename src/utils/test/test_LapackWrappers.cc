#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include "utils/LapackWrappers.h"
#include "utils/Utilities.h"

using namespace AMP;


// The main function
int main(int,char**) 
{
    int N_errors = 0;
    
    // Run the basic tests
    N_errors += Lapack::run_all_test();

    // Get the times for the tests
    const char* tests[] = { "dcopy", "dscal", "dnrm2", "dasum", "ddot", 
        "daxpy", "dgemv", "dgemm", "dgesv", "dgtsv", "dgbsv", 
        "dgetrf", "dgttrf","dgbtrf", "dgetrs", "dgttrs", "dgbtrs", 
        "dgetri" };
    const int N[] = { 500, 500, 500, 500, 500,
        100, 100, 100, 100, 500, 500,
        100, 500, 500, 100, 500, 500, 
        100 };
    for (size_t i=0; i<sizeof(tests)/sizeof(char*); i++) {
        double t1 = Utilities::time();
        int err = Lapack::run_test( tests[i], N[i] );
        double t2 = Utilities::time();
        int us = static_cast<int>(1e6*(t2-t1)/N[i]);
        printf("%7s:  %s:  %5i us\n",tests[i],err==0?"passed":"failed",us);
        N_errors += err;
    }

    // Finished
    if ( N_errors==0 )
        std::cout << "All tests passed\n";
    else
        std::cout << "Some tests failed\n";
    return N_errors==0 ? 0:1;
}


