#include "AMP/utils/AMPManager.h"
#include "AMP/AMP_Version.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/RNG.h"
#include "AMP/utils/ShutdownRegistry.h"
#include "AMP/utils/Utilities.h"

#include "LapackWrappers.h"
#include "ProfilerApp.h"
#include "StackTrace/ErrorHandlers.h"
#include "StackTrace/StackTrace.h"
#include "StackTrace/Utilities.h"


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
// clang-format on


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
#include <memory>


namespace AMP {


// Initialize static member variables
int AMPManager::initialized                 = 0;
int AMPManager::rank                        = 0;
bool AMPManager::called_PetscInitialize     = false;
bool AMPManager::use_MPI_Abort              = true;
bool AMPManager::print_times                = false;
int AMPManager::argc                        = 0;
char **AMPManager::argv                     = nullptr;
AMPManagerProperties AMPManager::properties = AMPManagerProperties();


/****************************************************************************
 *  Function to terminate AMP with a message for exceptions                  *
 ****************************************************************************/
static int force_exit      = 0;
static bool printed_stack  = false;
static int abort_stackType = 3;
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
    auto stack = StackTrace::getCallStack();
    for ( auto &elem : stack ) {
        if ( strcmp( elem.function.data(), "MPID_Abort" ) == 0 )
            return;
    }
    std::stringstream msg;
    msg << "Calling exit without calling shutdown\n";
    msg << "Bytes used = " << AMP::Utilities::getMemoryUsage() << std::endl;
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
        AMP_ERROR( "AMP was previously initialized.  It cannot be reinitialized" );
    if ( initialized == -1 )
        AMP_ERROR( "AMP was previously initialized and shutdown.  It cannot be reinitialized" );
    // Begin startup procedure
    double start    = Utilities::time();
    argc            = argc_in;
    argv            = argv_in;
    properties      = properties_in;
    print_times     = properties.print_times;
    abort_stackType = properties.stack_trace_type;
    // Initialize the timers (default is disabled)
    PROFILE_DISABLE();
    // Set the abort method
    AMPManager::use_MPI_Abort = properties.use_MPI_Abort;
    // Initialize MPI
    double MPI_start = Utilities::time();
    AMP_MPI::start_MPI( argc, argv, properties.profile_MPI_level );
    double MPI_time = Utilities::time() - MPI_start;
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
    StackTrace::Utilities::setErrorHandlers();
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
            printf( "  MPI startup time = %0.3f s\n", MPI_time );
        if ( petsc_time != 0 )
            printf( "  PETSc startup time = %0.3f s\n", petsc_time );
        if ( petsc_time != 0 )
            printf( "  SAMRAI startup time = %0.3f s\n", SAMRAI_time );
        printf( "\n" );
    }
}


/****************************************************************************
 * Shutdown the AMP package                                                  *
 ****************************************************************************/
void AMPManager::shutdown()
{
    double start_time = Utilities::time();
    int rank          = comm_world.getRank();
    if ( initialized == 0 )
        AMP_ERROR( "AMP is not initialized, did you forget to call startup or call shutdown more "
                   "than once" );
    if ( initialized == -1 )
        AMP_ERROR(
            "AMP has been initialized and shutdown.  Calling shutdown more than once is invalid" );
    initialized = -1;
    // Clear error handlers
    clearHandlers();
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
    double MPI_start = Utilities::time();
    AMP_MPI::stop_MPI();
    double MPI_time = Utilities::time() - MPI_start;
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
            printf( "  SAMRAI shutdown time = %0.3f s\n", SAMRAI_time );
        if ( petsc_time != 0 )
            printf( "  PETSc shutdown time = %0.3f s\n", petsc_time );
        if ( MPI_time != 0 )
            printf( "  MPI shutdown time = %0.3f s\n", MPI_time );
        printf( "\n" );
    }
    // Shutdown timer and print memory leaks on rank 0
    PROFILE_DISABLE();
#ifdef USE_TIMER
    auto memory = MemoryApp::getMemoryStats();
    if ( rank == 0 && memory.N_new > memory.N_delete )
        MemoryApp::print( std::cout );
#endif
    // Wait 50 milli-seconds for all processors to finish
    Utilities::sleepMs( 50 );
}


/****************************************************************************
 * Function to start/stop SAMRAI                                             *
 ****************************************************************************/
double AMPManager::start_SAMRAI()
{
    double time = 0;
#ifdef USE_EXT_SAMRAI
    double start = Utilities::time();
#ifdef USE_MPI
    SAMRAI::tbox::SAMRAI_MPI::init( AMP_MPI( AMP_COMM_WORLD ).getCommunicator() );
#else
    SAMRAI::tbox::SAMRAI_MPI::initMPIDisabled();
#endif
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
    SAMRAI::tbox::PIO::finalize();
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
class SAMRAIAbortAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    void logMessage( const std::string &msg, const std::string &file, const int line ) override
    {
        AMP::Utilities::abort( msg, file, line );
    }
    SAMRAIAbortAppender()          = default;
    virtual ~SAMRAIAbortAppender() = default;
};
#endif


/****************************************************************************
 * Function to start/stop PETSc                                              *
 ****************************************************************************/
double AMPManager::start_PETSc()
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
        PetscPopErrorHandler();
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
double AMPManager::start_CUDA()
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
    auto appender = std::make_shared<SAMRAIAbortAppender>();
    SAMRAI::tbox::Logger::getInstance()->setAbortAppender( appender );
#endif
    // Set the terminate routine for runtime errors
    StackTrace::Utilities::setErrorHandlers();
    // Set atexit function
    std::atexit( exitFun );
#ifdef USE_LINUX
    std::at_quick_exit( exitFun );
#endif
}
void AMPManager::clearHandlers()
{
    // Don't call the global version of the call stack
    StackTrace::globalCallStackFinalize();
    // Clear the MPI error handler for comm_world
    clearMPIErrorHandler();
    // Clear error handlers for StackTrace
    StackTrace::Utilities::clearErrorHandlers();
    StackTrace::clearSignals();
    StackTrace::clearSymbols();
    // Clear the error handlers for petsc
#ifdef USE_EXT_PETSC
    PetscPopSignalHandler();
    PetscPopErrorHandler();
#endif
}


/****************************************************************************
 *  Functions to handle MPI errors                                           *
 ****************************************************************************/
void AMPManager::setMPIErrorHandler()
{
#ifdef USE_EXT_MPI
    StackTrace::setMPIErrorHandler( comm_world.getCommunicator() );
#ifdef USE_EXT_SAMRAI
    StackTrace::setMPIErrorHandler( SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getCommunicator() );
#endif
#endif
}
void AMPManager::clearMPIErrorHandler()
{
#ifdef USE_EXT_MPI
    StackTrace::clearMPIErrorHandler( comm_world.getCommunicator() );
#ifdef USE_EXT_SAMRAI
    StackTrace::clearMPIErrorHandler(
        SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getCommunicator() );
#endif
#endif
}


/****************************************************************************
 * Empty constructor to setup default AMPManagerProperties                   *
 ****************************************************************************/
AMPManagerProperties::AMPManagerProperties()
    : use_MPI_Abort( true ),
      print_times( false ),
      profile_MPI_level( 2 ),
      print_startup( false ),
      stack_trace_type( 3 ),
      COMM_WORLD( AMP_COMM_WORLD )
{
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


/****************************************************************************
 *  Functions to return version info                                         *
 ****************************************************************************/
std::array<int, 3> AMPManager::revision()
{
    return { { AMP::Version::major, AMP::Version::minor, AMP::Version::build } };
}
std::string AMPManager::info()
{
    std::stringstream out;
    out << "AMP:" << std::endl;
    out << "   Version: " << AMP::Version::major << "." << AMP::Version::minor << "."
        << AMP::Version::build << std::endl;
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
    out << "SAMRAI: " << SAMRAI_VERSION_MAJOR << "." << SAMRAI_VERSION_MINOR << "."
        << SAMRAI_VERSION_PATCHLEVEL << std::endl;
#endif
#ifdef USE_EXT_PETSC
    out << "PETSc: " << PETSC_VERSION_MAJOR << "." << PETSC_VERSION_MINOR << "."
        << PETSC_VERSION_SUBMINOR << std::endl;
#endif
#ifdef USE_EXT_TRILINOS
    out << "Trilinos: " << TRILINOS_VERSION_STRING << std::endl;
#endif
#ifdef USE_EXT_SUNDIALS
    out << "Sundials: " << SUNDIALS_PACKAGE_VERSION << std::endl;
#endif
#ifdef HYPRE_RELEASE_VERSION
    out << "Hypre: " << HYPRE_RELEASE_VERSION << std::endl;
#elif defined( HYPRE_PACKAGE_VERSION )
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
    out << "MPI: " << AMP::AMP_MPI::info() << std::endl;
#endif
    out << "Lapack: " << Lapack<double>::info();
    return out.str();
}


} // namespace AMP
