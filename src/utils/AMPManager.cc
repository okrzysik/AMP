#include "AMPManager.h"
#include "PIO.h"
#include "ShutdownRegistry.h"
#include "Utilities.h"
#include "RNG.h"
#include "ProfilerApp.h"
#include "utils/AMP_MPI.h"

#ifdef USE_EXT_PETSC
    #include "petsc.h"
    #include "petscsys.h"
    #include "petscerror.h"
#endif
#ifdef USE_TIMER
    #include "MemoryApp.h"
#endif


#ifdef USE_EXT_LIBMESH
    //#include "ampmesh/libmesh/initializeLibMesh.h"
    /*#include "libmesh/enum_elem_quality.h"
    #include "libmesh/enum_elem_type.h"
    #include "libmesh/enum_eigen_solver_type.h"
    #include "libmesh/enum_fe_family.h"
    #include "libmesh/enum_inf_map_type.h"
    #include "libmesh/enum_io_package.h"
    #include "libmesh/enum_norm_type.h"
    #include "libmesh/enum_order.h"
    #include "libmesh/enum_parallel_type.h"
    #include "libmesh/enum_point_locator_type.h"
    #include "libmesh/enum_preconditioner_type.h"
    #include "libmesh/enum_quadrature_type.h"
    #include "libmesh/enum_solver_package.h"
    #include "libmesh/enum_solver_type.h"
    #include "libmesh/enum_subset_solve_mode.h"
    #include "libmesh/enum_xdr_mode.h"
    #include "libmesh/elem.h"
    #define define_libmesh_enum_map( TYPE, NAME )           \
        extern std::map<TYPE,std::string> enum_to_##NAME;   \
        extern std::map<std::string,TYPE> NAME##_to_enum;
    #define clear_libmesh_enum_map( TYPE, NAME )            \
        libMesh::enum_to_##NAME.clear();           \
        libMesh::NAME##_to_enum.clear();
    #include <map>
    #include <string>
    namespace libMesh {
        namespace {
            define_libmesh_enum_map(::ElemType,elem_type)
            define_libmesh_enum_map(::Order,order)
            define_libmesh_enum_map(::FEFamily,fefamily)
            define_libmesh_enum_map(::InfMapType,inf_map_type)
            define_libmesh_enum_map(::QuadratureType,quadrature_type)
            define_libmesh_enum_map(::PreconditionerType,preconditioner_type)
            define_libmesh_enum_map(::Elem::RefinementState,refinementstate_type)
            define_libmesh_enum_map(::EigenSolverType,eigensolvertype)
            define_libmesh_enum_map(::SolverType,solvertype)
            define_libmesh_enum_map(::ElemQuality,elemquality)
            define_libmesh_enum_map(::IOPackage,iopackage)
            define_libmesh_enum_map(::FEMNormType,norm_type)
            define_libmesh_enum_map(::ParallelType,parallel_type)
            define_libmesh_enum_map(::PointLocatorType,point_locator_type)
            define_libmesh_enum_map(::SolverPackage,solverpackage_type)
            define_libmesh_enum_map(::SubsetSolveMode,subset_solve_mode)
            define_libmesh_enum_map(::XdrMODE,xdr_mode)
        }
    }
    void clear_libmesh_enums() {
        clear_libmesh_enum_map(::ElemType,elem_type)
        clear_libmesh_enum_map(::Order,order)
        clear_libmesh_enum_map(::FEFamily,fefamily)
        clear_libmesh_enum_map(::InfMapType,inf_map_type)
        clear_libmesh_enum_map(::QuadratureType,quadrature_type)
        clear_libmesh_enum_map(::PreconditionerType,preconditioner_type)
        clear_libmesh_enum_map(::Elem::RefinementState,refinementstate_type)
        clear_libmesh_enum_map(::EigenSolverType,eigensolvertype)
        clear_libmesh_enum_map(::SolverType,solvertype)
        clear_libmesh_enum_map(::ElemQuality,elemquality)
        clear_libmesh_enum_map(::IOPackage,iopackage)
        clear_libmesh_enum_map(::FEMNormType,norm_type)
        clear_libmesh_enum_map(::ParallelType,parallel_type)
        clear_libmesh_enum_map(::PointLocatorType,point_locator_type)
        clear_libmesh_enum_map(::SolverPackage,solverpackage_type)
        clear_libmesh_enum_map(::SubsetSolveMode,subset_solve_mode)
        clear_libmesh_enum_map(::XdrMODE,xdr_mode)
    }*/
#endif

#include <new>
#include <string.h>
#include <iostream>
#include <sstream>

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
static int force_exit = 0;
void AMPManager::terminate_AMP(std::string message)
{
    if ( AMP::AMPManager::use_MPI_Abort==true || force_exit>0 ) {
        // Print the call stack and memory usage
        char text[100];
        std::stringstream msg;
        msg << message << std::endl;
        long long unsigned int N_bytes = AMP::Utilities::getMemoryUsage();
        sprintf(text,"Bytes used = %llu\n",N_bytes);
        msg << text;
        std::vector<std::string> stack = AMP::Utilities::getCallStack();
        msg << "Stack Trace:\n";
        for (size_t i=0; i<stack.size(); i++)
            msg << "   " << stack[i] << std::endl;
        perr << msg.str();
    }
    if ( force_exit>1 ) {
        exit(-1);
    } else if ( AMP::AMPManager::use_MPI_Abort==true ) {
        // Use MPI_abort (will terminate all processes)
        force_exit = 2;
        AMP_MPI(AMP_COMM_WORLD).abort();
    } else if ( force_exit>0 ) {
        exit(-1);
    } else {
        // Throw and standard exception (allows the use of try, catch)
        // std::stringstream  stream;
        // stream << message << std::endl << "  " << filename << ":  " << line;
        // std::cout << stream.str() << std::endl;
        throw std::logic_error(message);
    }
}


/****************************************************************************
*  Function to terminate AMP if an unhandled exception is caught            *
****************************************************************************/
void term_func_abort(int) 
{
    AMPManager::terminate_AMP("");
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
    force_exit = 1;
    AMPManager::terminate_AMP( "Unhandled exception:\n" + last_message );
}


/****************************************************************************
*  Function to handle MPI errors                                            *
****************************************************************************/
#ifdef USE_EXT_MPI
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
    AMPManager::terminate_AMP( msg.str() );
}
#endif


/****************************************************************************
*  Function to PETSc errors                                                 *
****************************************************************************/
#ifdef USE_EXT_PETSC
#if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==0 )
PetscErrorCode petsc_err_handler( int line, const char*, const char* file, 
    const char* dir, PetscErrorCode, int, const char* buf, void*)
#elif ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==2 )
PetscErrorCode petsc_err_handler( MPI_Comm, int line, const char*, 
    const char* dir, const char* file, PetscErrorCode, PetscErrorType, const char* buf, void* )
#else
#error Not programmed for this version of petsc
#endif
{
    std::stringstream msg;
    msg << "PETSc error:" << std::endl;
    msg << "   File: " << dir << file << ", line: " << line << std::endl;
    msg << "   " << buf << std::endl;
    AMPManager::terminate_AMP( msg.str());
    return 0;
}
#endif


/****************************************************************************
*                                                                            *
* Initialize the AMP package.  This routine performs the following tasks:   *
*                                                                            *
* (1) Initialize MPI                                                        *
*                                                                            *
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
            std::vector<char*> petscArgs = getPetscArgs();
            int n_args = static_cast<int>(petscArgs.size());
            char** args = NULL;
            if ( n_args>0 ) 
                args = &petscArgs[0];
            PetscInitialize(&n_args, &(args),  PETSC_NULL,  PETSC_NULL);
            called_PetscInitialize = true;
            for (size_t i=0; i<petscArgs.size(); i++)
                delete [] petscArgs[i];
        }
        // Set our error handler
        PetscPopSignalHandler();
        PetscPushErrorHandler(&petsc_err_handler,PETSC_NULL);
        petsc_time = time()-petsc_start_time;
    #endif
    // Initialize MPI
    AMP::AMP_MPI::changeProfileLevel( properties.profile_MPI_level );
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
    #endif
    setMPIErrorHandler();
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
*                                                                           *
* Shutdown the AMP package.  This routine currently only deallocates        *
* statically allocated memory and finalizes the output streams.             *
*                                                                           *
****************************************************************************/
void AMPManager::shutdown()
{    
    double start_time = time();
    double shutdown_time=0, MPI_time=0;
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
    #ifdef USE_EXT_LIBMESH
        //if ( AMP::Mesh::initializeLibMesh::isInitialized() ) {
        //    AMP_ERROR("Libmesh should be finalized before shutting down");
        //}
        // delete std::map objects from libmesh string_to_enum.C
        //clear_libmesh_enums();
    #endif
    // Shutdown MPI
    comm_world.barrier();
    AMPManager::use_MPI_Abort = false;
    comm_world = AMP_MPI(AMP_COMM_NULL);    // Delete comm world
    clearMPIErrorHandler();
    if ( called_MPI_Init ) {
        double MPI_start_time = time();
        #ifdef USE_EXT_MPI
            MPI_Finalize();
        #endif
        MPI_time = time()-MPI_start_time;
    }
    // Shudown PETSc
    double petsc_time = 0.0;
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
    // Print any AMP_MPI leaks
    if ( AMP_MPI::MPI_Comm_created()!=AMP_MPI::MPI_Comm_destroyed() ) {
        printf("Rank %i detected AMP_MPI comm leak: %i %i\n",rank,
            static_cast<int>(AMP_MPI::MPI_Comm_created()),
            static_cast<int>(AMP_MPI::MPI_Comm_destroyed()));
    }
    // Shutdown the profiler
    PROFILE_DISABLE();
    // Print memory leaks on rank 0
    #ifdef USE_TIMER
        MemoryApp::MemoryStats memory = MemoryApp::getMemoryStats();
        if ( rank==0 && memory.N_new>memory.N_delete )
            MemoryApp::print(std::cout);
    #endif
    // Wait 50 milli-seconds for all processors to finish
    Sleep(50);
}


/****************************************************************************
* Function to create the arguments for petsc                                *
****************************************************************************/
static inline void addArg( std::string arg, std::vector<char*>& args)
{
    char *tmp = new char[arg.length()+1];
    memset( tmp, 0, arg.length()+1 );
    memcpy( tmp, arg.c_str(), arg.length() );
    args.push_back(tmp);
}
std::vector<char*> AMPManager::getPetscArgs()
{
    std::vector<char*> args;
    addArg( "-malloc no", args);
    return args;
}


/****************************************************************************
* Functions to set/clear the MPI error handler                              *
****************************************************************************/
#ifdef USE_EXT_MPI
    boost::shared_ptr<MPI_Errhandler> AMPManager::mpierr;
#endif
void AMPManager::setMPIErrorHandler( )
{
    #ifdef USE_EXT_MPI
        if ( mpierr.get()==NULL ) {
            mpierr = boost::shared_ptr<MPI_Errhandler>( new MPI_Errhandler );
            MPI_Comm_create_errhandler( MPI_error_handler_fun, mpierr.get() );
        }
        MPI_Comm_set_errhandler( MPI_COMM_SELF, *mpierr );
        MPI_Comm_set_errhandler( MPI_COMM_WORLD, *mpierr );
    #endif
}
void AMPManager::clearMPIErrorHandler(  )
{
    #ifdef USE_EXT_MPI
        if ( mpierr.get()!=NULL )
            MPI_Errhandler_free( mpierr.get() );    // Delete the error handler
        mpierr.reset();
        MPI_Comm_set_errhandler( MPI_COMM_SELF, MPI_ERRORS_ARE_FATAL );
        MPI_Comm_set_errhandler( MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL );
    #endif
}


/****************************************************************************
* Empty constructor to setup default AMPManagerProperties                   *
****************************************************************************/
AMPManagerProperties::AMPManagerProperties() 
{
    use_MPI_Abort = true;
    print_times = false;
    profile_MPI_level = 2;
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

