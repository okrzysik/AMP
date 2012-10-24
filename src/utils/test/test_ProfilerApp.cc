#include "utils/AMPManager.h"
#include "utils/ProfilerApp.h"




int main(int argc, char* argv[])
{
    // Initialize AMP
    AMP::AMPManager::startup(argc,argv);
    PROFILE_ENABLE();

    PROFILE_START("MAIN");
    int k;
    for (int i=0; i<1e5; i++) {
        PROFILE_START("gettimeofday");
        timeval time1;
        k = 0;
        for (int j=0; j<100; j++) {
            gettimeofday(&time1,NULL);
            k++;
        }
        PROFILE_STOP("gettimeofday");
        /*PROFILE_START("gethrtime");
        hrtime_t time2;
        k = 0;
        for (int j=0; j<10; j++) {
            time2 = gethrtime(); 
            k++;
        }
        PROFILE_STOP("gethrtime");*/
    }
    PROFILE_STOP("MAIN");
    PROFILE_SAVE("test_ProfilerApp.timer");

    // Finalize AMP
    AMP::AMPManager::shutdown();
    return(0);
}



