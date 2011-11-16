#include "AMPManager.h"
#include "PIO.h"
#include "ShutdownRegistry.h"
#include "Utilities.h"
#include "RNG.h"
//#include "materials/Material.h"

#ifdef USE_MPI
    #include "mpi.h"
#endif

#ifdef USE_PETSC
    #include "petscsys.h"   
#endif

#ifdef USE_LIBMESH
    #include "ampmesh/libmesh/initializeLibMesh.h"
#endif

#include <new>
#include <string.h>


#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    // Windows 
    // Sleep is defined in milliseconds
#else
    // Linux
    // usleep is defined in microseconds, create a Sleep command
    #define Sleep(x) usleep(x*1000)
#endif


namespace AMP {


// Initialize static member variables
bool AMPManager::initialized=false;
bool AMPManager::called_MPI_Init=false;
bool AMPManager::called_PetscInitialize=false;
bool AMPManager::use_MPI_Abort=true;
bool AMPManager::print_times=false;
AMP_MPI AMPManager::comm=AMP::AMP_MPI();
int AMPManager::argc=0;
char** AMPManager::argv=NULL;
AMPManagerProperties AMPManager::properties=AMPManagerProperties();


// Function to get the current time (preferably using a hi resolution timer
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    #include <windows.h>
    double time() { 
        LARGE_INTEGER end, f;
        QueryPerformanceFrequency(&f);
        QueryPerformanceCounter(&end);       
        double time = ((double)end.QuadPart)/((double)f.QuadPart));
        return time;
    }
#else
    #include <sys/time.h>
    double time() { 
        timeval current_time;
        gettimeofday(&current_time,NULL);
        double time = ((double)current_time.tv_sec)+1e-6*((double)current_time.tv_usec);
        return time;
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
    if ( initialized )
        AMP_ERROR("AMP was previously initialized");
    double start_time = time();
    double startup_time=0, petsc_time=0, MPI_time=0;
    argc = argc_in;
    argv = argv_in;
    properties = properties_in;
    print_times = properties.print_times;
    // Set the abort method
    AMPManager::use_MPI_Abort = properties.use_MPI_Abort;
    // Initialize PETSc
    #ifdef USE_PETSC
        double petsc_start_time = time();
        if ( PetscInitializeCalled ) {
            called_PetscInitialize = false;
        } else {
            PetscInitialize(&argc, &argv,  PETSC_NULL,  PETSC_NULL);
            called_PetscInitialize = true;
        }
        petsc_time = time()-petsc_start_time;
    #endif
    // Initialize MPI
    #ifdef USE_MPI
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
    #endif
    // Initialize AMP's MPI
    comm = AMP_MPI(AMP_COMM_WORLD);
    int size1 = comm.getSize();
    #ifdef USE_MPI
        int size2 = 0;
        MPI_Comm_size(MPI_COMM_WORLD,&size2);
    #else
        int size2 = 1;
    #endif
    AMP_INSIST(size1==size2,"AMP's MPI size does not match MPI_COMM_WORLD");
    // Initialize the parallel IO
    PIO::initialize();
    // Initialize the random number generator
    AMP::RNG::initialize(123);
    // Initialize the Materials interface
    //AMP::Materials::initialize();
    // Initialization finished
    initialized = true;
    startup_time = time()-start_time;
    if ( print_times && comm.getRank()==0 ) {
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
    int rank = comm.getRank();
    // Check if AMP was previously initialized
    if ( !initialized )
        AMP_ERROR("AMP is not initialized, did you forget to call startup or call shutdown more than once");
    // Syncronize all processors
    comm.barrier();
    ShutdownRegistry::callRegisteredShutdowns();
    // Shutdown the parallel IO
    PIO::finalize();
    // Shutdown LibMesh
    #ifdef USE_LIBMESH
        if ( AMP::Mesh::initializeLibMesh::isInitialized() ) {
            AMP_ERROR("Libmesh should be finalized before shutting down");
        }
    #endif
    // Shutdown MPI
    if ( called_MPI_Init ) {
        double MPI_start_time = time();
        MPI_Finalize();
        MPI_time = time()-MPI_start_time;
    }
    // Shudown PETSc
    #ifdef USE_PETSC
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

