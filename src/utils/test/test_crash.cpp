// This tests the error handler making sure we catch segfaults and print the call stack
#include "AMP/IO/PIO.h"
#include "AMP/utils/AMPManager.h"

#include <csignal>
#include <iomanip>
#include <iostream>
#include <sstream>


void crash()
{
    AMP::AMP_MPI comm( AMP_COMM_WORLD );
    if ( comm.getRank() == 0 )
        raise( SIGSEGV );
    comm.barrier();
}


int main( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::logOnlyNodeZero( "test_crash.log" );
    crash();
    printf( "FAILED" );
    AMP::AMPManager::shutdown();
    return 1;
}
