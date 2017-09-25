// This tests the error handler making sure we catch segfaults and print the call stack
#include "utils/AMPManager.h"
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
    crash();
    printf( "FAILED" );
    AMP::AMPManager::shutdown();
    return 1;
}
