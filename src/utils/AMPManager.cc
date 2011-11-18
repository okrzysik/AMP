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
    #include "libmesh.h"
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
AMP_MPI AMPManager::comm_world=AMP::AMP_MPI();
void* AMPManager::lminit=NULL;
int AMPManager::argc_libmesh=0;
char** AMPManager::argv_libmesh=NULL;


// Function to alter the command line arguments for libmesh
#ifdef USE_LIBMESH
static int add_libmesh_cmdline ( int argc , char **argv, char ***argv_new )
{
    const int N_add = 0;    // Number of additional arguments we want to add
    // Copy the existing command-line arguments (shifting by the number of additional arguments)
    *argv_new = new char*[argc+N_add];
    for ( int i = 0 ; i != argc ; i++ ) {
        (*argv_new)[i] = new char [ strlen ( argv[i] ) + 1 ];
        strcpy ( (*argv_new)[i] , argv[i] );
    }
    /*// Add command to keep cout from all processors (not just rank 0)
    (*argv_new)[argc] = new char [12];
    strcpy ( (*argv_new)[argc] , "--keep-cout" );*/
    return argc+N_add;
}
#endif


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
void AMPManager::startup(int argc, char *argv[], const AMPManagerProperties &properties)
{
    double start_time = time();
    double startup_time=0, petsc_time=0, MPI_time=0, libmesh_time=0;
    print_times = properties.print_times;
    // Check if AMP was previously initialized
    if ( initialized )
        AMP_ERROR("AMP was previously initialized");
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
    if ( properties.COMM_WORLD == AMP_COMM_WORLD ) 
        comm_world = AMP_MPI(MPI_COMM_WORLD);
    else
        comm_world = AMP_MPI(properties.COMM_WORLD);
    // Initialize LibMesh (must be done after initializing PETSc and MPI)
    #ifdef USE_LIBMESH
        double libmesh_start_time = time();
        argc_libmesh = add_libmesh_cmdline(argc,argv,&argv_libmesh);
        lminit = (void*) new LibMeshInit(argc_libmesh,argv_libmesh);
        // Reset the PETSc communicator back to the global communicator
        // LibMesh resets it to it's communicator regardless of whether it initializes it or not
        #ifdef USE_PETSC
            PETSC_COMM_WORLD = MPI_COMM_WORLD;
        #endif
        libmesh_time = time()-libmesh_start_time;
    #endif
    // Initialize the parallel IO
    PIO::initialize();
    // Initialize the random number generator
    AMP::RNG::initialize(123);
    // Initialize the Materials interface
    //AMP::Materials::initialize();
    // Initialization finished
    initialized = true;
    startup_time = time()-start_time;
    if ( print_times && comm_world.getRank()==0 ) {
        printf("startup time = %0.3f s\n",startup_time);
        if ( petsc_time!=0 )
            printf(" PETSc startup time = %0.3f s\n",petsc_time);
        if ( libmesh_time!=0 )
            printf(" libmesh startup time = %0.3f s\n",libmesh_time);
        if ( MPI_time!=0 )
            printf(" libmesh startup time = %0.3f s\n",MPI_time);
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
    double shutdown_time=0, petsc_time=0, MPI_time=0, libmesh_time=0;
    int rank = comm_world.getRank();
    // Check if AMP was previously initialized
    if ( !initialized )
        AMP_ERROR("AMP is not initialized, did you forget to call startup or call shutdown more than once");
    // Syncronize all processors
    comm_world.barrier();
    ShutdownRegistry::callRegisteredShutdowns();
    // Shutdown the parallel IO
    PIO::finalize();
    // Shtudown LibMesh
    #ifdef USE_LIBMESH
        double libmesh_start_time = time();
        delete (LibMeshInit*) lminit;
        libmesh_time = time()-libmesh_start_time;
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
    // Delete internal variables
    if ( argv_libmesh!=NULL ) {
        for (int i=0; i<argc_libmesh; i++)
            delete [] argv_libmesh[i];
        delete[] argv_libmesh;
        argv_libmesh = NULL;
    }
    argc_libmesh = 0;
    Sleep(10);
    shutdown_time = time()-start_time;
    if ( print_times && rank==0 ) {
        printf("shutdown time = %0.3f s\n",shutdown_time);
        if ( petsc_time!=0 )
            printf(" PETSc shutdown time = %0.3f s\n",petsc_time);
        if ( MPI_time!=0 )
            printf(" libmesh shutdown time = %0.3f s\n",MPI_time);
        if ( libmesh_time!=0 )
            printf(" libmesh shutdown time = %0.3f s\n",libmesh_time);
    }
    // Wait 50 milli-seconds for all processors to finish
    initialized = false;
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
*									                                        *
* Function to initialize Libmesh with a given communicator                  *
* We want to preinitialize LibMesh with the global communicator used to     *
* initialize PETSc and MPI.   Unfortunatelly, LibMesh needs to be           *
* reinitialized with a different communicator if we are using a load        *
* balancing that distributes the meshes into sub-groups of processors.      *
* Since LibMesh does not have a clean way to do this and connot have        *
* multiple LibMeshInit this function is provided to re-initialize libmesh.  *
*									                                        *
****************************************************************************/
void AMPManager::initializeLibmesh( AMP_MPI libmeshComm )
{
    #ifdef USE_LIBMESH
        // Delete the existing LibMeshInit object so we can re-initialize LibMesh
        // Note:  This only works if MPI and PETSc have been initialized outside of LibMesh
        delete (LibMeshInit*) lminit;
        // Reinitialize LibMesh with the new communicator
        #ifdef USE_MPI
            AMP::AMPManager::lminit = new LibMeshInit ( argc_libmesh, argv_libmesh, libmeshComm.getCommunicator() );
        #else
            AMP::AMPManager::lminit = new LibMeshInit ( argc_libmesh, argv_libmesh );
        #endif
        // Reset the PETSc communicator back to the global communicator
        // LibMesh resets it to it's communicator regardless of weather it initializes it or not
        #ifdef USE_MPI
            PETSC_COMM_WORLD = MPI_COMM_WORLD;
        #endif
    #else
        AMP_ERROR("Libmesh is disables, initializeLibmesh(AMP_MPI comm) is non-functional");
    #endif
}



}

