#include "AMP/utils/AMPManager.h"
#include "AMP/AMP_Version.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/RNG.h"
#include "AMP/utils/ShutdownRegistry.h"
#include "AMP/utils/Utilities.h"

#include "LapackWrappers.h"
#include "ProfilerApp.h"
#include "StackTrace/StackTrace.h"

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
#ifdef USE_EXT_SAMRAI
    #undef NULL_USE
    #include "SAMRAI/tbox/StartupShutdownManager.h"
    #include "SAMRAI/tbox/SAMRAIManager.h"
    #include "SAMRAI/tbox/Logger.h"
#endif
#ifdef USE_EXT_HYPRE
    #include "HYPRE_config.h"
#endif
#ifdef USE_CUDA
    #include <cuda_runtime_api.h>
    #include <cuda.h>
#endif
// clang-format off


#include <iostream>
#include <new>
#include <sstream>
#include <cstring>
#include <vector>
#include <array>
#include <algorithm>
#include <csignal>
#include <stdexcept>
#include <cstdio>


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
static int abort_stackType = 3;
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
        StackTrace::multi_stack_info stack;
        if ( abort_stackType == 1 ) {
            stack = StackTrace::getCallStack();
        } else if ( abort_stackType == 2 ) {
            stack = StackTrace::getAllCallStacks();
        } else if ( abort_stackType == 3 ) {
            stack = StackTrace::getGlobalCallStacks();
        }
        StackTrace::cleanupStackTrace( stack );
        auto data = stack.print();
        msg << std::endl;
        msg << "Stack Trace:\n";
        for ( const auto &i : data )
            msg << " " << i << std::endl;
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
    auto stack = StackTrace::getCallStack();
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
* Initialize the AMP package.  This routine performs the following tasks:   *
****************************************************************************/
void AMPManager::startup( int argc_in, char *argv_in[], const AMPManagerProperties &properties_in )
{
    // Check if AMP was previously initialized
    if ( initialized == 1 )
        AMP_ERROR( "AMP was previously initialized and shutdown.  It cannot be reinitialized" );
    if ( initialized == -1 )
        AMP_ERROR( "AMP was previously initialized and shutdown.  It cannot be reinitialized" );
    // Begin startup procedure
    double start = Utilities::time();
    argc         = argc_in;
    argv         = argv_in;
    properties   = properties_in;
    print_times  = properties.print_times;
    // Initialize the timers (default is disabled)
    PROFILE_DISABLE();
    // Set the abort method
    AMPManager::use_MPI_Abort = properties.use_MPI_Abort;
    // Initialize MPI
    double MPI_time = start_MPI( argc, argv, properties.profile_MPI_level );
    // Initialize AMP's MPI
    if ( properties.COMM_WORLD == AMP_COMM_WORLD )
        comm_world = AMP_MPI( MPI_COMM_WORLD );
    else
        comm_world = AMP_MPI( properties.COMM_WORLD );
    // Initialize PETSc
    double petsc_time = start_PETSc();
    // Initialize SAMRAI
    double SAMRAI_time = start_SAMRAI();
    // Initialize the parallel IO
    PIO::initialize();
    // Initialze call stack
    StackTrace::globalCallStackInitialize( comm_world.getCommunicator() );
    // Initialize the random number generator
    AMP::RNG::initialize( 123 );
    // Initialize cuda
    start_CUDA();
    // Set the signal/terminate handlers
    setHandlers();
    // Initialization finished
    initialized = 1;
    rank        = comm_world.getRank();
    double time = Utilities::time() - start;
    if ( properties.print_startup ) {
        AMP::pout << "Version info:\n" << info() << std::endl;
        AMP::pout.flush();
    }
    if ( print_times && comm_world.getRank() == 0 ) {
        printf( "startup time = %0.3f s\n", time );
        if ( MPI_time != 0 )
            printf( " MPI startup time = %0.3f s\n", MPI_time );
        if ( petsc_time != 0 )
            printf( " PETSc startup time = %0.3f s\n", petsc_time );
        if ( petsc_time != 0 )
            printf( " SAMRAI startup time = %0.3f s\n", SAMRAI_time );
        printf( "\n" );
    }
}


/****************************************************************************
* Shutdown the AMP package                                                  *
****************************************************************************/
void AMPManager::shutdown()
{
    double start_time = Utilities::time();
    int rank = comm_world.getRank();
    if ( initialized == 0 )
        AMP_ERROR( "AMP is not initialized, did you forget to call startup or call shutdown more "
                   "than once" );
    if ( initialized == -1 )
        AMP_ERROR(
            "AMP has been initialized and shutdown.  Calling shutdown more than once is invalid" );
    initialized = -1;
    // Disable call stack and error handlers
    StackTrace::globalCallStackFinalize( );
    clearMPIErrorHandler();
    // Disable MPI_Abort
    AMPManager::use_MPI_Abort = false;
    // Shutdown the registry
    ShutdownRegistry::callRegisteredShutdowns();
    // Shutdown the parallel IO
    PIO::finalize();
    // Shutdown the profiler
    PROFILE_DISABLE();
    // Syncronize all ranks
    comm_world.barrier();
    // Shudown PETSc
    double petsc_time = stop_PETSc();
    // Shutdown SAMRAI
    double SAMRAI_time = stop_SAMRAI();
    // Shutdown MPI
    double MPI_time = stop_MPI();
    // Print any AMP_MPI leaks
    if ( AMP_MPI::MPI_Comm_created() != AMP_MPI::MPI_Comm_destroyed() ) {
        printf( "Rank %i detected AMP_MPI comm leak: %i %i\n",
                rank,
                static_cast<int>( AMP_MPI::MPI_Comm_created() ),
                static_cast<int>( AMP_MPI::MPI_Comm_destroyed() ) );
    }
    // List shutdown times
    double shutdown_time = Utilities::time() - start_time;
    if ( print_times && rank == 0 ) {
        printf( "shutdown time = %0.3f s\n", shutdown_time );
        if ( SAMRAI_time != 0 )
            printf( " SAMRAI shutdown time = %0.3f s\n", SAMRAI_time );
        if ( petsc_time != 0 )
            printf( " PETSc shutdown time = %0.3f s\n", petsc_time );
        if ( MPI_time != 0 )
            printf( " MPI shutdown time = %0.3f s\n", MPI_time );
        printf( "\n" );
    }
// Print memory leaks on rank 0
#ifdef USE_TIMER
    MemoryApp::MemoryStats memory = MemoryApp::getMemoryStats();
    if ( rank == 0 && memory.N_new > memory.N_delete )
        MemoryApp::print( std::cout );
#endif
    // Wait 50 milli-seconds for all processors to finish
    Utilities::sleepMs( 50 );
}


/****************************************************************************
* Function to start/stop MPI                                                *
****************************************************************************/
double AMPManager::start_MPI( int argc, char *argv[], int profile_level )
{
    double time = 0;
    AMP::AMP_MPI::changeProfileLevel( profile_level );
    NULL_USE( argc );
    NULL_USE( argv );
#ifdef USE_EXT_MPI
    if ( MPI_Active() ) {
        called_MPI_Init = false;
    } else {
        int provided;
        double start = Utilities::time();
        int result   = MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided );
        if ( result != MPI_SUCCESS )
            AMP_ERROR( "AMP was unable to initialize MPI" );
        if ( provided < MPI_THREAD_MULTIPLE )
            AMP::perr << "Warning: Failed to start MPI with MPI_THREAD_MULTIPLE\n";
        called_MPI_Init = true;
        time            = Utilities::time() - start;
    }
#endif
    return time;
}
double AMPManager::stop_MPI()
{
    double time = 0;
    comm_world = AMP_MPI( AMP_COMM_NULL );
#ifdef USE_EXT_MPI
    int finalized;
    MPI_Finalized( &finalized );
    if ( called_MPI_Init && !finalized ) {
        double start = Utilities::time();
        MPI_Finalize();
        time = Utilities::time() - start;
    }
#endif
    return time;
}


/****************************************************************************
* Function to start/stop SAMRAI                                             *
****************************************************************************/
double AMPManager::start_SAMRAI()
{
    double time = 0;
#ifdef USE_EXT_SAMRAI
    double start = Utilities::time();
    SAMRAI::tbox::SAMRAI_MPI::init( AMP_MPI( AMP_COMM_WORLD ).getCommunicator() );
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();
    SAMRAI::tbox::SAMRAIManager::setMaxNumberPatchDataEntries( 2048 );
    time = Utilities::time() - start;
#endif
    return time;
}
double AMPManager::stop_SAMRAI()
{
    double time = 0;
#ifdef USE_EXT_SAMRAI
    double start = Utilities::time();
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();
    time = Utilities::time() - start;
#endif
    return time;
}


/****************************************************************************
 *  Class to override the output appender for abort messages                 *
 ****************************************************************************/
#ifdef USE_EXT_SAMRAI
#include "SAMRAI/tbox/Logger.h"
class SAMRUtilsAbortAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    void
    logMessage( const std::string &msg, const std::string &file, const int line ) override
    {
        AMP::Utilities::abort( msg, file, line );
    }
    SAMRUtilsAbortAppender() = default;
    virtual ~SAMRUtilsAbortAppender() = default;
};
#endif


/****************************************************************************
* Function to start/stop PETSc                                              *
****************************************************************************/
double AMPManager::start_PETSc( )
{
    double time = 0;
#ifdef USE_EXT_PETSC
    double start = Utilities::time();
    if ( PetscInitializeCalled ) {
        called_PetscInitialize = false;
    } else {
        auto petscArgs = getPetscArgs();
        auto n_args    = static_cast<int>( petscArgs.size() );
        char **args    = nullptr;
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
    time = Utilities::time() - start;
#endif
    return time;
}
double AMPManager::stop_PETSc()
{
    double time = 0;
#ifdef USE_EXT_PETSC
    if ( called_PetscInitialize ) {
        double start = Utilities::time();
        PetscPopSignalHandler();
        PetscFinalize();
        time = Utilities::time() - start;
    }
#endif
    return time;
}


/****************************************************************************
* Function to create the arguments for petsc                                *
****************************************************************************/
static inline void addArg( std::string arg, std::vector<char *> &args )
{
    auto *tmp = new char[arg.length() + 1];
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
* Function to start/stop CUDA                                               *
****************************************************************************/
double AMPManager::start_CUDA( )
{
#ifdef USE_CUDA
    void *tmp;
    cudaMallocManaged( &tmp, 10, cudaMemAttachGlobal );
    cudaFree( tmp );
#endif
    return 0;
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
    PetscPopErrorHandler();
    PetscPushErrorHandler( &petsc_err_handler, PETSC_NULL );
#endif
// Set the error handlers for SAMRAI
#ifdef USE_EXT_SAMRAI
    SAMRAI::tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit( true );
    SAMRAI::tbox::SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort( true );
    auto appender = AMP::make_shared<SAMRUtilsAbortAppender>();
    SAMRAI::tbox::Logger::getInstance()->setAbortAppender( appender );
#endif
    // Set the terminate routine for runtime errors
    StackTrace::setErrorHandlers( abort_fun );
    // Set atexit function
    std::atexit( exitFun );
#ifdef USE_LINUX
    std::at_quick_exit( exitFun );
#endif
}
void dummyErrorHandler( std::string, StackTrace::terminateType ) {}
void AMPManager::clearHandlers()
{
    initialized = 3;
    // Don't call the global version of the call stack
    StackTrace::globalCallStackFinalize();
    // Clear the MPI error handler for comm_world
    clearMPIErrorHandler();
    // Clear the error handlers for petsc
#ifdef USE_EXT_PETSC
    PetscPopSignalHandler();
    PetscPopErrorHandler();
#endif
}


/****************************************************************************
*  Functions to handle MPI errors                                           *
****************************************************************************/
#ifdef USE_EXT_MPI
AMP::shared_ptr<MPI_Errhandler> AMPManager::mpierr;
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
void AMPManager::setMPIErrorHandler( )
{
#ifdef USE_EXT_MPI
    if ( MPI_Active() ) {
        if ( mpierr.get() == nullptr ) {
            mpierr = AMP::make_shared<MPI_Errhandler>( );
            MPI_Comm_create_errhandler( MPI_error_handler_fun, mpierr.get() );
        }
        MPI_Comm_set_errhandler( MPI_COMM_SELF, *mpierr );
        MPI_Comm_set_errhandler( MPI_COMM_WORLD, *mpierr );
        if ( comm_world.getCommunicator() != MPI_COMM_WORLD )
            MPI_Comm_set_errhandler( comm_world.getCommunicator(), *mpierr );
#ifdef USE_EXT_SAMRAI
        MPI_Comm_set_errhandler( SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getCommunicator(), *mpierr );
#endif
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
#ifdef USE_EXT_SAMRAI
        MPI_Comm_set_errhandler( SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getCommunicator(), MPI_ERRORS_ARE_FATAL );
#endif
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
    return { {AMP::Version::major, AMP::Version::minor, AMP::Version::build} };
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
#ifdef USE_EXT_SAMRAI
    out << "SAMRAI: " << SAMRAI_VERSION_MAJOR << "." <<
                        SAMRAI_VERSION_MINOR << "." <<
                        SAMRAI_VERSION_PATCHLEVEL << std::endl;
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
#ifdef HYPRE_RELEASE_VERSION
    out << "Hypre: " << HYPRE_RELEASE_VERSION << std::endl;
#elif defined(HYPRE_PACKAGE_VERSION)
    out << "Hypre: " << HYPRE_PACKAGE_VERSION << std::endl;
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
        out << "MPI: " << MPI_info << std::endl;
    } else {
        int MPI_version;
        int MPI_subversion;
        MPI_Get_version( &MPI_version, &MPI_subversion );
        out << "MPI: " << MPI_version << "." << MPI_subversion << std::endl;
    }
#endif
    out << "Lapack: " << Lapack<double>::info();
    return out.str();
}


} // namespace AMP
