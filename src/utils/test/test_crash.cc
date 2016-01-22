// This tests the error handler making sure we catch segfaults and print the call stack
#include "utils/AMPManager.h"
#include <iomanip>
#include <iostream>
#include <signal.h>
#include <sstream>


void crash() { raise( SIGSEGV ); }


int main( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    crash();
    printf( "FAILED" );
    AMP::AMPManager::shutdown();
    return 1;
}
