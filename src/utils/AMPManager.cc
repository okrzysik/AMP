#include "AMPManager.h"
#include "PIO.h"
#include "ShutdownRegistry.h"
#include "Utilities.h"
#include "RNG.h"
#include "utils/ProfilerApp.h"
#include "utils/AMP_MPI.h"

#ifdef USE_EXT_PETSC
    #include "petsc.h"
    #include "petscsys.h"
    #include "petscerror.h"
#endif

//#ifdef USE_EXT_LIBMESH
//    #include "ampmesh/libmesh/initializeLibMesh.h"
//#endif

#include <new>
#include <string.h>
#include <iostream>
#include <sstream>


#ifdef _MSC_VER
    #define _CRT_SECURE_NO_WARNINGS		// Supress depreciated warnings for visual studio
#endif

#include <stdio.h>
#include <stdexcept>
#include <signal.h>


namespace AMP {


// Initialize static member variables
int AMPManager::initialized=0;
bool AMPManager::called_MPI_Init=false;
bool AMPManager::called_PetscInitialize=false;
bool AMPManager::use_MPI_Abort=true;
bool AMPManager::print_times=false;
AMP_MPI AMPManager::comm_world=AMP::AMP_MPI();
int AMPManager::argc=0;
char** AMPManager::argv=NULL;
AMPManagerProperties AMPManager::properties=AMPManagerProperties();


// Function to get the current time (preferably using a hi resolution timer
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    #define USE_WINDOWS
    #include <windows.h>
    // Sleep is defined in milliseconds
    double time() { 
        LARGE_INTEGER end, f;
        QueryPerformanceFrequency(&f);
        QueryPerformanceCounter(&end);       
        double time = ((double)end.QuadPart)/((double)f.QuadPart);
        return time;
    }
#else
    #define USE_LINUX
    #include <sys/time.h>
    // usleep is defined in microseconds, create a Sleep command
    #define Sleep(x) usleep(x*1000)
    double time() { 
        timeval current_time;
        gettimeofday(&current_time,NULL);
        double time = ((double)current_time.tv_sec)+1e-6*((double)current_time.tv_usec);
        return time;
    }
#endif


/****************************************************************************
*  Function to terminate AMP with a message for exceptions                  *
****************************************************************************/
void terminate_AMP(std::string message)
{
    std::stringstream msg;
    char text[100];
    msg << message << std::endl;
    long long unsigned int N_bytes = AMP::Utilities::getMemoryUsage();
    sprintf(text,"Bytes used = %llu\n",N_bytes);
    msg << text;
    std::vector<std::string> stack = AMP::Utilities::getCallStack();
    msg << "Stack Trace:\n";
    for (size_t i=0; i<stack.size(); i++)
        msg << "   " << stack[i];
    std::cerr << msg.str();
    AMP_MPI(AMP_COMM_WORLD).abort();
}


/****************************************************************************
*  Function to terminate AMP if an unhandled exception is caught            *
****************************************************************************/
void term_func_abort(int err) 
{
    terminate_AMP("");
}
static int tried_throw = 0;
void term_func() 
{
    // Try to re-throw the last error to get the last message
    std::string last_message;
    #ifdef USE_LINUX
        try {
            if ( tried_throw==0 ) { 
                tried_throw = 1;
                throw;
            }
            // No active exception
        } catch (const std::exception &err) {
            // Caught a std::runtime_error
            last_message = err.what();
        } catch (...) {
            // Caught an unknown exception
            last_message = "unknown exception occurred.";
        }
    #endif
    terminate_AMP( "Unhandled exception:\n" + last_message );
}


/****************************************************************************
*  Function to handle MPI errors                                            *
****************************************************************************/
#ifdef USE_EXT_MPI
MPI_Errhandler AMPManager::mpierr;
static void MPI_error_handler_fun( MPI_Comm *comm, int *err, ... )
{
    if ( *err==MPI_ERR_COMM && *comm==MPI_COMM_WORLD ) {
        // Special error handling for an invalid MPI_COMM_WORLD
        std::cerr << "Error invalid MPI_COMM_WORLD";
        exit(-1);
    }
    int msg_len=0;
    std::stringstream msg;    
    char message[1000];
    MPI_Error_string( *err, message, &msg_len );
    if ( msg_len <= 0 )
        AMP_ERROR("Unkown error in MPI");
    msg << "Error calling MPI routine:\n" + std::string(message) << std::endl;
    terminate_AMP( msg.str() );
}
#endif


/****************************************************************************
*  Function to PETSc errors                                                 *
****************************************************************************/
#ifdef USE_EXT_PETSC
PetscErrorCode petsc_err_handler(int line, const char* func, const char* file, const char* dir, PetscErrorCode code, int p, const char* buf, void* ctx)
{
    std::stringstream msg;
    msg << "PETSc error:" << std::endl;
    msg << "   File: " << dir << file << ", line: " << line << std::endl;
    msg << "   " << buf << std::endl;
    terminate_AMP( msg.str());
    return 0;
}
#endif


/****************************************************************************
*									                                        *
* Initialize the AMP package.  This routine performs the following tasks:   *
*									                                        *
* (1) Initialize MPI                                                        *
*									                                        *
****************************************************************************/
void AMPManager::startup(int argc_in, char *argv_in[], const AMPManagerProperties &properties_in)
{
    // Check if AMP was previously initialized
    if ( initialized==1 )
        AMP_ERROR("AMP was previously initialized and shutdown.  It cannot be reinitialized");
    if ( initialized==-1 )
        AMP_ERROR("AMP was previously initialized and shutdown.  It cannot be reinitialized");
    double start_time = time();
    double startup_time=0, petsc_time=0, MPI_time=0;
    argc = argc_in;
    argv = argv_in;
    properties = properties_in;
    print_times = properties.print_times;
    // Initialize the timers (default is disabled)
    PROFILE_DISABLE();
    // Set the abort method
    AMPManager::use_MPI_Abort = properties.use_MPI_Abort;
    // Initialize PETSc
    #ifdef USE_EXT_PETSC
        double petsc_start_time = time();
        if ( PetscInitializeCalled ) {
            called_PetscInitialize = false;
        } else {
            PetscInitialize(&argc, &argv,  PETSC_NULL,  PETSC_NULL);
            called_PetscInitialize = true;
        }
        // Set our error handler
        PetscPopSignalHandler();
        PetscPushErrorHandler(&petsc_err_handler,PETSC_NULL);
        petsc_time = time()-petsc_start_time;
    #endif
    // Initialize MPI
    #ifdef USE_EXT_MPI
        int flag;
        MPI_Initialized(&flag);
        if ( flag ) {
            called_MPI_Init = false;
            MPI_time = 0;
        } else {
            double MPI_start_time = time();
            int result = MPI_Init(&argc, &argv);
            if (result != MPI_SUCCESS) 
                AMP_ERROR("AMP was unable to initialize MPI");
            called_MPI_Init = true;
            MPI_time = time()-MPI_start_time;
        }
        MPI_Comm_create_errhandler( MPI_error_handler_fun, &mpierr );
    #endif
    // Initialize AMP's MPI
    if ( properties.COMM_WORLD == AMP_COMM_WORLD ) 
		#ifdef USE_EXT_MPI
			comm_world = AMP_MPI(MPI_COMM_WORLD);
		#else
			comm_world = AMP_MPI(AMP_COMM_WORLD);
		#endif
    else
        comm_world = AMP_MPI(properties.COMM_WORLD);    // Initialize the parallel IO
    PIO::initialize();
    // Initialize the random number generator
    AMP::RNG::initialize(123);
    // Set the terminate routine for runtime errors
    std::set_terminate( term_func );
    signal(SIGABRT,&term_func_abort);
    signal(SIGSEGV,&term_func_abort);
    //std::set_unexpected( term_func );
    // Initialization finished
    initialized = 1;
    startup_time = time()-start_time;
    if ( print_times && comm_world.getRank()==0 ) {
        printf("startup time = %0.3f s\n",startup_time);
        if ( petsc_time!=0 )
            printf(" PETSc startup time = %0.3f s\n",petsc_time);
         if ( MPI_time!=0 )
            printf(" MPI startup time = %0.3f s\n",MPI_time);
    }
}


/****************************************************************************
*									                                        *
* Shutdown the AMP package.  This routine currently only deallocates	    *
* statically allocated memory and finalizes the output streams.		        *
*									                                        *
****************************************************************************/
void AMPManager::shutdown()
{    
    double start_time = time();
    double shutdown_time=0, petsc_time=0, MPI_time=0;
    int rank = comm_world.getRank();
    if ( initialized==0 )
        AMP_ERROR("AMP is not initialized, did you forget to call startup or call shutdown more than once");
    if ( initialized==-1 )
        AMP_ERROR("AMP has been initialized and shutdown.  Calling shutdown more than once is invalid");
    // Syncronize all processors
    comm_world.barrier();
    ShutdownRegistry::callRegisteredShutdowns();
    // Shutdown the parallel IO
    PIO::finalize();
    // Shutdown LibMesh
    /*#ifdef USE_EXT_LIBMESH
        if ( AMP::Mesh::initializeLibMesh::isInitialized() ) {
            AMP_ERROR("Libmesh should be finalized before shutting down");
        }
    #endif*/
    // Shutdown MPI
    comm_world = AMP_MPI(AMP_COMM_NULL);    // Delete comm world
    #ifdef USE_EXT_MPI
        MPI_Errhandler_free( &mpierr );    // Delete the error handler
        MPI_Comm_set_errhandler( MPI_COMM_SELF, MPI_ERRORS_ARE_FATAL );
        MPI_Comm_set_errhandler( MPI_COMM_SELF, MPI_ERRORS_ARE_FATAL );
    #endif
    if ( called_MPI_Init ) {
        double MPI_start_time = time();
        #ifdef USE_EXT_MPI
            MPI_Finalize();
        #endif
        MPI_time = time()-MPI_start_time;
    }
    // Shudown PETSc
    #ifdef USE_EXT_PETSC
        if ( called_PetscInitialize ) {
            double petsc_start_time = time();
            PetscFinalize();
            petsc_time = time()-petsc_start_time;
        }
    #endif
    Sleep(10);
    shutdown_time = time()-start_time;
    if ( print_times && rank==0 ) {
        printf("shutdown time = %0.3f s\n",shutdown_time);
        if ( petsc_time!=0 )
            printf(" PETSc shutdown time = %0.3f s\n",petsc_time);
        if ( MPI_time!=0 )
            printf(" MPI shutdown time = %0.3f s\n",MPI_time);
    }
    // Wait 50 milli-seconds for all processors to finish
    Sleep(50);
}



/****************************************************************************
*									                                        *
* Empty constructor to setup default AMPManagerProperties                   *
*									                                        *
****************************************************************************/
AMPManagerProperties::AMPManagerProperties() {
    use_MPI_Abort = true;
    print_times = false;
    COMM_WORLD = AMP_COMM_WORLD;
}


/****************************************************************************
*  Some simple functions                                                    *
****************************************************************************/
int AMPManager::get_argc() {
    AMP_INSIST(initialized,"AMP has not been initialized");
    return argc;
}
char** AMPManager::get_argv() {
    AMP_INSIST(initialized,"AMP has not been initialized");
    return argv;
}
AMPManagerProperties AMPManager::getAMPManagerProperties() {
    AMP_INSIST(initialized,"AMP has not been initialized");
    return properties;
}



}

