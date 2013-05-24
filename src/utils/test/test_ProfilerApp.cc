#include "utils/AMPManager.h"
#include "utils/ProfilerApp.h"

#include <string>

#ifdef USE_WINDOWS
    #define TIME_TYPE LARGE_INTEGER
    #define get_time(x) QueryPerformanceCounter(x)
    #define get_diff(start,end,f) (((double)(end.QuadPart-start.QuadPart))/((double)f.QuadPart))
    #define get_frequency(f) QueryPerformanceFrequency(f)
#elif defined(USE_LINUX) || defined(USE_MAC)
    #define TIME_TYPE timeval
    #define get_time(x) gettimeofday(x,NULL);
    #define get_diff(start,end,f) (((double)end.tv_sec-start.tv_sec)+1e-6*((double)end.tv_usec-start.tv_usec))
    #define get_frequency(f) (*f=timeval())
#else
    #error Unknown OS
#endif


int run_tests( bool enable_trace, std::string save_name ) 
{
    PROFILE_ENABLE();
    PROFILE_SYNCRONIZE();
    if ( enable_trace ) {
        PROFILE_ENABLE_TRACE();
        PROFILE_ENABLE_MEMORY();
    }
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
        // Test how long it takes to get the time
        PROFILE_START("gettime");
        TIME_TYPE time1;
        for (int j=0; j<N_timers; j++)
            get_time(&time1);
        PROFILE_STOP("gettime");
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
        // Test the memory around allocations
        PROFILE_START("allocate1");
        PROFILE_START("allocate2");
        double *tmp = new double[5000000];
        NULL_USE(tmp);
        PROFILE_STOP("allocate2");
        delete [] tmp;
        PROFILE_START("allocate3");
        tmp = new double[100000];
        NULL_USE(tmp);
        PROFILE_STOP("allocate3");
        delete [] tmp;
        PROFILE_STOP("allocate1");
    }

    PROFILE_STOP("MAIN");
    PROFILE_SAVE(save_name);
    return N_errors;
}


int main(int argc, char* argv[])
{
    // Initialize AMP
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);
    
    // Run the tests
    int N_errors=0;
    N_errors += run_tests( false, "test_ProfilerApp" );
    PROFILE_DISABLE();
    N_errors += run_tests( true, "test_ProfilerApp-trace" );
    
    // Finalize AMP
    if ( N_errors==0 ) 
        std::cout << "All tests passed" << std::endl;
    else
        std::cout << "Some tests failed" << std::endl;
    AMP::AMPManager::shutdown();
    return(N_errors);
}



