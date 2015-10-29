#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include "utils/LapackWrappers.h"
#include "utils/Utilities.h"

using namespace AMP;


// The main function
int main( int, char*[] ) 
{
    int N_errors = 0;
    
    // Print the machine specifics
    Lapack::print_machine_parameters( );

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
        double error = 0;
        int err = Lapack::run_test( tests[i], N[i], error );
        double t2 = Utilities::time();
        int us = static_cast<int>(1e6*(t2-t1)/N[i]);
        printf("%7s:  %s:  %5i us  (%e)\n",tests[i],err==0?"passed":"failed",us,error);
        N_errors += err;
    }

    // Run the tests in parallel to check for parallel bugs
    printf("\nRunning parallel tests\n");
    int N_threads = 8;
    ThreadPool tpool(N_threads);
    std::vector<double> error(N_threads,0);
    std::vector<ThreadPool::WorkItem*> work(N_threads,NULL);
    std::vector<ThreadPool::thread_id_t> id(N_threads);
    for (size_t i=0; i<sizeof(tests)/sizeof(char*); i++) {
        double t1 = Utilities::time();
        for (int j=0; j<N_threads; j++)
            work[j] = new WorkItemFull<int,const char*,int,double&>(
                Lapack::run_test, tests[i], N[i], error[j] );
        tpool.add_work(work);
        tpool.wait_pool_finished();
        double t2 = Utilities::time();
        bool pass = true;
        for (int j=0; j<N_threads; j++)
            pass = pass && tpool.getFunctionRet<int>(id[j])==0;
        if ( !pass ) {
            printf("Failed parallel test: %s\n",tests[i]);
            N_errors++;
        }
        int us = static_cast<int>(1e6*(t2-t1)/(N[i]*N_threads));
        printf("%7s:  %5i us\n",tests[i],us);
    }


    // Finished
    if ( N_errors==0 )
        std::cout << "\nAll tests passed\n";
    else
        std::cout << "\nSome tests failed\n";
    return N_errors==0 ? 0:1;
}


