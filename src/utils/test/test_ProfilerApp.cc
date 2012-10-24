#include "utils/AMPManager.h"
#include "utils/ProfilerApp.h"

#include <string>


inline void test_timer( const std::string& name, int level ) {
    PROFILE_START(name,level);
    PROFILE_STOP(name,level);
}



int main(int argc, char* argv[])
{
    // Initialize AMP
    AMP::AMPManager::startup(argc,argv);
    PROFILE_ENABLE();

    const int N_it = 100;
    const int N_timers = 1000;

    // Get a list of timer names
    std::vector<std::string> names(N_timers);
    for (int i=0; i<N_timers; i++) {
        char tmp[16];
        sprintf(tmp,"%04i",i);
        names[i] = std::string(tmp);
    }

    PROFILE_START("MAIN");
    int k;
    for (int i=0; i<1e2; i++) {
        // Test how long it takes to get the time of day
        PROFILE_START("gettimeofday");
        timeval time1;
        k = 0;
        for (int j=0; j<N_timers; j++) {
            gettimeofday(&time1,NULL);
            k++;
        }
        PROFILE_STOP("gettimeofday");
        // Test how long it takes to start/stop the timers
        PROFILE_START("level 0");
        for (int j=0; j<N_timers; j++)
            test_timer(names[j],0);
        PROFILE_STOP("level 0");
        PROFILE_START("level 1");
        for (int j=0; j<N_timers; j++)
            test_timer(names[j],1);
        PROFILE_STOP("level 1");
    }
    PROFILE_STOP("MAIN");
    PROFILE_SAVE("test_ProfilerApp");

    // Finalize AMP
    AMP::AMPManager::shutdown();
    return(0);
}



