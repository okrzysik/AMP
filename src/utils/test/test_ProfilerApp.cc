#include "utils/AMPManager.h"
#include "utils/ProfilerApp.h"

#include <string>



int main(int argc, char* argv[])
{
    // Initialize AMP
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);
    PROFILE_ENABLE();
    PROFILE_START("MAIN");

    const int N_it = 100;
    const int N_timers = 1000;
    int N_errors = 0;

    // Get a list of timer names
    std::vector<std::string> names(N_timers);
    for (int i=0; i<N_timers; i++) {
        char tmp[16];
        sprintf(tmp,"%04i",i);
        names[i] = std::string(tmp);
    }

    // Check that the start/stop command fail when they should
    try {   // Check basic start/stop
        PROFILE_START("dummy1"); 
        PROFILE_STOP("dummy1"); 
    } catch (... ) {
        N_errors++;
    }
    try {   // Check stop call before start
        PROFILE_STOP("dummy2"); 
        N_errors++;
    } catch (... ) {
    }
    try {   // Check multiple calls to start
        PROFILE_START("dummy3"); 
        PROFILE_START2("dummy3"); 
        N_errors++;
    } catch (... ) {
        PROFILE_STOP("dummy3"); 
    }
    try {   // Check multiple calls to start with different line numbers
        PROFILE_START("dummy1");
        N_errors++;
    } catch (... ) {
    }

    // Check the performance
    for (int i=0; i<N_it; i++) {
        // Test how long it takes to get the time of day
        PROFILE_START("gettimeofday");
        timeval time1;
        for (int j=0; j<N_timers; j++)
            gettimeofday(&time1,NULL);
        PROFILE_STOP("gettimeofday");
        // Test how long it takes to start/stop the timers
        PROFILE_START("level 0");
        for (int j=0; j<N_timers; j++) {
            PROFILE_START(names[j],0);
            PROFILE_STOP(names[j],0);
        }
        PROFILE_STOP("level 0");
        PROFILE_START("level 1");
        for (int j=0; j<N_timers; j++) {
            PROFILE_START(names[j],1);
            PROFILE_STOP(names[j],1);
        }
        PROFILE_STOP("level 1");
    }


    PROFILE_STOP("MAIN");
    PROFILE_SAVE("test_ProfilerApp");

    // Finalize AMP
    if ( N_errors==0 ) 
        std::cout << "All tests passed" << std::endl;
    else
        std::cout << "Some tests failed" << std::endl;
    AMP::AMPManager::shutdown();
    return(N_errors);
}



