#include "AMPManager.h"
#include "AMP_Version.h"
#include "PIO.h"
#include "ProfilerApp.h"
#include "RNG.h"
#include "ShutdownRegistry.h"
#include "utils/AMP_MPI.h"
#include "utils/StackTrace.h"
#include "utils/Utilities.h"


// Include external files for startup/version info
// clang-format off
#ifdef USE_EXT_PETSC
    #include "petsc.h"
    #include "petscerror.h"
    #include "petscsys.h"
    #include "petscversion.h"
#endif
#if USE_EXT_TIMER
    #include "TimerUtilityVersion.h"
    #include "MemoryApp.h"
#endif
#ifdef USE_EXT_TRILINOS
    #include "Trilinos_version.h"
#endif
#ifdef USE_EXT_LIBMESH
    #include "libmesh/libmesh_version.h"
#endif
#ifdef USE_EXT_HDF5
    #include "H5public.h"
#endif
#ifdef USE_EXT_SUNDIALS
    #include "sundials/sundials_config.h"
#endif
#ifdef USE_EXT_SILO
    #include "silo.h"
#endif
#ifdef USE_EXT_HYPRE
    #include "HYPRE_config.h"
#endif
// clang-format off


#include <iostream>
#include <new>
#include <sstream>
#include <string.h>

#include <algorithm>
#include <signal.h>
#include <stdexcept>
#include <stdio.h>


namespace AMP {


// Initialize static member variables
int AMPManager::initialized                 = 0;
int AMPManager::rank                        = 0;
bool AMPManager::called_MPI_Init            = false;
bool AMPManager::called_PetscInitialize     = false;
bool AMPManager::use_MPI_Abort              = true;
bool AMPManager::print_times                = false;
AMP_MPI AMPManager::comm_world              = AMP::AMP_MPI();
int AMPManager::argc                        = 0;
char **AMPManager::argv                     = nullptr;
AMPManagerProperties AMPManager::properties = AMPManagerProperties();


/****************************************************************************
*  Function to terminate AMP with a message for exceptions                  *
****************************************************************************/
static int force_exit     = 0;
static bool printed_stack = false;
static void abort_fun( std::string msg, StackTrace::terminateType type )
{
    if ( type == StackTrace::terminateType::exception )
        force_exit = std::max( force_exit, 1 );
    AMPManager::terminate_AMP( msg );
}
void AMPManager::terminate_AMP( std::string message )
{
    AMP_MPI comm( AMP_COMM_WORLD );
    if ( !printed_stack ) {
        // Print the call stack and memory usage
        std::stringstream msg;
        msg << message << std::endl;
        msg << "Bytes used = " << AMP::Utilities::getMemoryUsage() << std::endl;
        auto stack = AMP::StackTrace::getCallStack();
        msg << "Stack Trace:\n";
        for ( auto &elem : stack )
            msg << "   " << elem.print() << std::endl;
        // Add a rank dependent wait to hopefully print the stack trace cleanly
        Utilities::sleepMs( ( 100 * comm.getRank() ) / comm.getSize() );
        perr << msg.str();
        printed_stack = true;
        force_exit    = 1;
    }
    if ( force_exit > 1 ) {
        exit( -1 );
    } else if ( AMP::AMPManager::use_MPI_Abort == true ) {
        // Use MPI_abort (will terminate all processes)
        force_exit = 2;
        comm.abort();
    } else if ( force_exit > 0 ) {
        exit( -1 );
    } else {
        // Throw and standard exception (allows the use of try, catch)
        force_exit = 1;
        throw std::logic_error( message );
    }
}
void AMPManager::exitFun()
{
    if ( initialized != 1 || printed_stack )
        return;
    std::stringstream msg;
    msg << "Calling exit without calling shutdown\n";
    msg << "Bytes used = " << AMP::Utilities::getMemoryUsage() << std::endl;
    auto stack = AMP::StackTrace::getCallStack();
    msg << "Stack Trace:\n";
    for ( auto &elem : stack )
        msg << "   " << elem.print() << std::endl;
    perr << msg.str();
}


/****************************************************************************
*  Function to PETSc errors                                                 *
****************************************************************************/
#ifdef USE_EXT_PETSC
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
PetscErrorCode petsc_err_handler( int line,
                                  const char *,
                                  const char *file,
                                  const char *dir,
                                  PetscErrorCode,
                                  int,
                                  const char *buf,
                                  void * )
#elif ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 2 )
PetscErrorCode petsc_err_handler( MPI_Comm,
                                  int line,
                                  const char *,
                                  const char *dir,
                                  const char *file,
                                  PetscErrorCode,
                                  PetscErrorType,
                                  const char *buf,
                                  void * )
#elif ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 7 )
PetscErrorCode petsc_err_handler( MPI_Comm,
                                  int line,
                                  const char *dir,
                                  const char *file,
                                  PetscErrorCode,
                                  PetscErrorType,
                                  const char *buf,
                                  void * )
#else
#error Not programmed for this version of petsc
#endif
{
    std::stringstream msg;
    msg << "PETSc error:" << std::endl;
    msg << "   File: " << dir << file << ", line: " << line << std::endl;
    msg << "   " << buf << std::endl;
    AMPManager::terminate_AMP( msg.str() );
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
void AMPManager::startup( int argc_in, char *argv_in[], const AMPManagerProperties &properties_in )
{
    // Check if AMP was previously initialized
    if ( initialized == 1 )
        AMP_ERROR( "AMP was previously initialized and shutdown.  It cannot be reinitialized" );
    if ( initialized == -1 )
        AMP_ERROR( "AMP was previously initialized and shutdown.  It cannot be reinitialized" );
    double start_time   = Utilities::time();
    double startup_time = 0, petsc_time = 0, MPI_time = 0;
    argc        = argc_in;
    argv        = argv_in;
    properties  = properties_in;
    print_times = properties.print_times;
    // Initialize the timers (default is disabled)
    PROFILE_DISABLE();
    // Set the abort method
    AMPManager::use_MPI_Abort = properties.use_MPI_Abort;
    // Initialize PETSc
#ifdef USE_EXT_PETSC
    double petsc_start_time = Utilities::time();
    if ( PetscInitializeCalled ) {
        called_PetscInitialize = false;
    } else {
        std::vector<char *> petscArgs = getPetscArgs();
        int n_args                    = static_cast<int>( petscArgs.size() );
        char **args                   = nullptr;
        if ( n_args > 0 )
            args = &petscArgs[0];
        PetscInitialize( &n_args, &( args ), PETSC_NULL, PETSC_NULL );
        called_PetscInitialize = true;
        for ( auto &petscArg : petscArgs )
            delete[] petscArg;
    }
    #ifndef USE_MPI
        // Fix minor bug in petsc where first call to dup returns MPI_COMM_WORLD instead of a new comm
        AMP::AMP_MPI( MPI_COMM_WORLD ).dup();
    #endif
    petsc_time = Utilities::time() - petsc_start_time;
#endif
    // Initialize MPI
    AMP::AMP_MPI::changeProfileLevel( properties.profile_MPI_level );
#ifdef USE_EXT_MPI
    if ( MPI_Active() ) {
        called_MPI_Init = false;
        MPI_time        = 0;
    } else {
        double MPI_start_time = Utilities::time();
        int result            = MPI_Init( &argc, &argv );
        if ( result != MPI_SUCCESS )
            AMP_ERROR( "AMP was unable to initialize MPI" );
        called_MPI_Init = true;
        MPI_time        = Utilities::time() - MPI_start_time;
    }
#endif
    // Initialize AMP's MPI
    if ( properties.COMM_WORLD == AMP_COMM_WORLD )
        comm_world = AMP_MPI( MPI_COMM_WORLD );
    else
        comm_world = AMP_MPI( properties.COMM_WORLD ); // Initialize the parallel IO
    PIO::initialize();
    // Initialize the random number generator
    AMP::RNG::initialize( 123 );
    // Set the signal/terminate handlers
    setHandlers();
    // Initialization finished
    initialized  = 1;
    rank         = comm_world.getRank();
    startup_time = Utilities::time() - start_time;
    if ( properties.print_startup ) {
        AMP::pout << "Version info:\n" << info() << std::endl;
        AMP::pout.flush();
    }
    if ( print_times && comm_world.getRank() == 0 ) {
        printf( "startup time = %0.3f s\n", startup_time );
        if ( petsc_time != 0 )
            printf( " PETSc startup time = %0.3f s\n", petsc_time );
        if ( MPI_time != 0 )
            printf( " MPI startup time = %0.3f s\n", MPI_time );
        printf( "\n" );
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
    double start_time    = Utilities::time();
    double shutdown_time = 0, MPI_time = 0;
    int rank = comm_world.getRank();
    if ( initialized == 0 )
        AMP_ERROR( "AMP is not initialized, did you forget to call startup or call shutdown more "
                   "than once" );
    if ( initialized == -1 )
        AMP_ERROR(
            "AMP has been initialized and shutdown.  Calling shutdown more than once is invalid" );
    initialized = -1;
    // Syncronize all processors
    comm_world.barrier();
    ShutdownRegistry::callRegisteredShutdowns();
    // Shutdown the parallel IO
    PIO::finalize();
// Shutdown LibMesh
#ifdef USE_EXT_LIBMESH
// if ( AMP::Mesh::initializeLibMesh::isInitialized() ) {
//    AMP_ERROR("Libmesh should be finalized before shutting down");
//}
// delete std::map objects from libmesh string_to_enum.C
// clear_libmesh_enums();
#endif
    // Shutdown MPI
    comm_world.barrier();   // Sync all processes
    clearMPIErrorHandler(); // Clear MPI error handler before deleting comms
    AMPManager::use_MPI_Abort = false;
    comm_world                = AMP_MPI( AMP_COMM_NULL ); // Delete comm world
    if ( called_MPI_Init ) {
        double MPI_start_time = Utilities::time();
#ifdef USE_EXT_MPI
        MPI_Finalize();
#endif
        MPI_time = Utilities::time() - MPI_start_time;
    }
    // Shudown PETSc
    double petsc_time = 0.0;
#ifdef USE_EXT_PETSC
    if ( called_PetscInitialize ) {
        double petsc_start_time = Utilities::time();
        PetscFinalize();
        petsc_time = Utilities::time() - petsc_start_time;
    }
#endif
    Utilities::sleepMs( 10 );
    shutdown_time = Utilities::time() - start_time;
    if ( print_times && rank == 0 ) {
        printf( "shutdown time = %0.3f s\n", shutdown_time );
        if ( petsc_time != 0 )
            printf( " PETSc shutdown time = %0.3f s\n", petsc_time );
        if ( MPI_time != 0 )
            printf( " MPI shutdown time = %0.3f s\n", MPI_time );
        printf( "\n" );
    }
    // Print any AMP_MPI leaks
    if ( AMP_MPI::MPI_Comm_created() != AMP_MPI::MPI_Comm_destroyed() ) {
        printf( "Rank %i detected AMP_MPI comm leak: %i %i\n",
                rank,
                static_cast<int>( AMP_MPI::MPI_Comm_created() ),
                static_cast<int>( AMP_MPI::MPI_Comm_destroyed() ) );
    }
    // Shutdown the profiler
    PROFILE_DISABLE();
// Print memory leaks on rank 0
#ifdef USE_TIMER
    MemoryApp::MemoryStats memory = MemoryApp::getMemoryStats();
    if ( rank == 0 && memory.N_new > memory.N_delete )
        MemoryApp::print( std::cout );
#endif
// Wait 50 milli-seconds for all processors to finish
#ifdef USE_EXT_MPI
    if ( MPI_Active() )
        MPI_Barrier( MPI_COMM_WORLD );
#endif
    Utilities::sleepMs( 50 );
}


/****************************************************************************
* Function to create the arguments for petsc                                *
****************************************************************************/
static inline void addArg( std::string arg, std::vector<char *> &args )
{
    char *tmp = new char[arg.length() + 1];
    memset( tmp, 0, arg.length() + 1 );
    memcpy( tmp, arg.c_str(), arg.length() );
    args.push_back( tmp );
}
std::vector<char *> AMPManager::getPetscArgs()
{
    std::vector<char *> args;
    addArg( "-malloc no", args );
    return args;
}


/****************************************************************************
* Functions to set/clear the error handlers                                 *
****************************************************************************/
void AMPManager::setHandlers()
{
    // Set the MPI error handler for comm_world
    setMPIErrorHandler();
// Set the error handlers for petsc
#ifdef USE_EXT_PETSC
    PetscPopSignalHandler();
    PetscPushErrorHandler( &petsc_err_handler, PETSC_NULL );
#endif
    // Set the terminate routine for runtime errors
    StackTrace::setErrorHandlers( abort_fun );
    // Set atexit function
    std::atexit( exitFun );
#ifdef USE_LINUX
    std::at_quick_exit( exitFun );
#endif
}
#ifdef USE_EXT_MPI
AMP::shared_ptr<MPI_Errhandler> AMPManager::mpierr;
#endif


/****************************************************************************
*  Functions to handle MPI errors                                           *
****************************************************************************/
#ifdef USE_EXT_MPI
static void MPI_error_handler_fun( MPI_Comm *comm, int *err, ... )
{
    if ( *err == MPI_ERR_COMM && *comm == MPI_COMM_WORLD ) {
        // Special error handling for an invalid MPI_COMM_WORLD
        std::cerr << "Error invalid MPI_COMM_WORLD";
        exit( -1 );
    }
    int msg_len = 0;
    char message[1000];
    MPI_Error_string( *err, message, &msg_len );
    if ( msg_len <= 0 )
        AMP_ERROR( "Unkown error in MPI" );
    std::string msg = "Error calling MPI routine:\n" + std::string( message ) + "\n";
    throw std::logic_error( msg );
}
#endif
void AMPManager::setMPIErrorHandler()
{
#ifdef USE_EXT_MPI
    if ( MPI_Active() ) {
        if ( mpierr.get() == nullptr ) {
            mpierr = AMP::shared_ptr<MPI_Errhandler>( new MPI_Errhandler );
            MPI_Comm_create_errhandler( MPI_error_handler_fun, mpierr.get() );
        }
        MPI_Comm_set_errhandler( MPI_COMM_SELF, *mpierr );
        MPI_Comm_set_errhandler( MPI_COMM_WORLD, *mpierr );
        if ( comm_world.getCommunicator() != MPI_COMM_WORLD )
            MPI_Comm_set_errhandler( comm_world.getCommunicator(), *mpierr );
    }
#endif
}
void AMPManager::clearMPIErrorHandler()
{
#ifdef USE_EXT_MPI
    if ( MPI_Active() ) {
        if ( mpierr.get() != nullptr )
            MPI_Errhandler_free( mpierr.get() ); // Delete the error handler
        mpierr.reset();
        MPI_Comm_set_errhandler( MPI_COMM_SELF, MPI_ERRORS_ARE_FATAL );
        MPI_Comm_set_errhandler( MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL );
        if ( comm_world.getCommunicator() != MPI_COMM_WORLD )
            MPI_Comm_set_errhandler( comm_world.getCommunicator(), MPI_ERRORS_ARE_FATAL );
    }
#endif
}


/****************************************************************************
* Empty constructor to setup default AMPManagerProperties                   *
****************************************************************************/
AMPManagerProperties::AMPManagerProperties()
{
    use_MPI_Abort     = true;
    print_times       = false;
    profile_MPI_level = 2;
    print_startup     = false;
    COMM_WORLD        = AMP_COMM_WORLD;
}


/****************************************************************************
*  Some simple functions                                                    *
****************************************************************************/
int AMPManager::get_argc()
{
    AMP_INSIST( initialized, "AMP has not been initialized" );
    return argc;
}
char **AMPManager::get_argv()
{
    AMP_INSIST( initialized, "AMP has not been initialized" );
    return argv;
}
AMPManagerProperties AMPManager::getAMPManagerProperties()
{
    AMP_INSIST( initialized, "AMP has not been initialized" );
    return properties;
}
bool AMPManager::MPI_Active()
{
#ifdef USE_EXT_MPI
    int MPI_initialized, MPI_finialized;
    MPI_Initialized( &MPI_initialized );
    MPI_Finalized( &MPI_finialized );
    return MPI_initialized != 0 && MPI_finialized == 0;
#else
    return false;
#endif
}


/****************************************************************************
*  Functions to return version info                                         *
****************************************************************************/
std::array<int,3> AMPManager::revision()
{
    return { AMP::Version::major, AMP::Version::minor, AMP::Version::build };
}
std::string AMPManager::info()
{
    std::stringstream out; 
    out << "AMP:" << std::endl;
    out << "   Version: " << AMP::Version::major << "." <<
                              AMP::Version::minor << "." <<
                              AMP::Version::build << std::endl;
    out << "   Hash: " << AMP::Version::short_hash << std::endl;
    out << "   C Compiler: " << AMP::Version::C << std::endl;
    out << "   C++ Compiler: " << AMP::Version::CXX << std::endl;
    out << "   Fortran Compiler: " << AMP::Version::Fortran << std::endl;
    out << "   C Compiler ID: " << AMP::Version::C_ID << std::endl;
    out << "   C++ Compiler ID: " << AMP::Version::CXX_ID << std::endl;
    out << "   Fortran Compiler ID: " << AMP::Version::Fortran_ID << std::endl;
    out << "   C Compiler Version: " << AMP::Version::C_VERSION << std::endl;
    out << "   C++ Compiler Version: " << AMP::Version::CXX_VERSION << std::endl;
    out << "   Fortran Compiler Version: " << AMP::Version::Fortran_VERSION << std::endl;
    out << "   C Flags: " << AMP::Version::C_FLAGS << std::endl;
    out << "   C++ Flags: " << AMP::Version::CXX_FLAGS << std::endl;
    out << "   Fortran Flags: " << AMP::Version::Fortran_FLAGS << std::endl;
#ifdef USE_TIMER
    out << "ProfilerApp: " << TIMER_VERSION << std::endl;
#endif
#ifdef USE_EXT_PETSC
    out << "PETSc: " << PETSC_VERSION_MAJOR << "." <<
                        PETSC_VERSION_MINOR << "." <<
                        PETSC_VERSION_SUBMINOR << std::endl;
#endif
#ifdef USE_EXT_TRILINOS
    out << "Trilinos: " << TRILINOS_VERSION_STRING << std::endl;
#endif
#ifdef USE_EXT_SUNDIALS
    out << "Sundials: " << SUNDIALS_PACKAGE_VERSION << std::endl;
#endif
#ifdef USE_EXT_HYPRE
    out << "Hypre: " << HYPRE_RELEASE_VERSION << std::endl;
#endif
#ifdef USE_EXT_LIBMESH
    out << "libMesh: " << libMesh::get_libmesh_version() << std::endl;
#endif
#ifdef USE_EXT_HDF5
    out << "HDF5: " << H5_VERS_MAJOR << "." << H5_VERS_MINOR << "." << H5_VERS_RELEASE << std::endl;
#endif
#ifdef USE_EXT_SILO
    out << "SILO: " << SILO_VERS_MAJ << "." << SILO_VERS_MIN << "." << SILO_VERS_PAT << std::endl;
#endif
#ifdef USE_EXT_MPI
    int MPI_version_length = 0;
    char MPI_version_string[MPI_MAX_LIBRARY_VERSION_STRING];
#if MPI_VERSION >= 3
    MPI_Get_library_version( MPI_version_string, &MPI_version_length );
#endif
    if ( MPI_version_length > 0 ) {
        std::string MPI_info( MPI_version_string, MPI_version_length );
        size_t pos = MPI_info.find('\n');
        while ( pos != std::string::npos ) {
             MPI_info.insert( pos+1, "   " );
             pos = MPI_info.find('\n',pos+1);
        }
        out << "MPI: " << MPI_info;
    } else {
        int MPI_version;
        int MPI_subversion;
        MPI_Get_version( &MPI_version, &MPI_subversion );
        out << "MPI: " << MPI_version << "." << MPI_subversion << std::endl;
    }
#endif
    return out.str();
}


} // namespace AMP
