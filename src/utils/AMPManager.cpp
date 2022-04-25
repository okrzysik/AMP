#include "AMP/utils/AMPManager.h"
#include "AMP/AMP_TPLs.h"
#include "AMP/AMP_Version.h"
#include "AMP/IO/PIO.h"
#include "AMP/utils/AMP_MPI.I"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/RNG.h"
#include "AMP/utils/Utilities.h"

#include "ProfilerApp.h"
#include "StackTrace/ErrorHandlers.h"
#include "StackTrace/StackTrace.h"
#include "StackTrace/Utilities.h"


// Include external files for startup/version info
// clang-format off
#ifdef USE_CUDA
    #include <cuda.h>
    #include <cuda_runtime_api.h>
#endif
#ifdef AMP_USE_LAPACK_WRAPPERS
    #include "LapackWrappers.h"
#endif
#ifdef AMP_USE_PETSC
    #include "petsc.h"
    #include "petscerror.h"
    #include "petscsys.h"
    #include "petscversion.h"
#endif
#ifdef AMP_USE_TIMER
    #include "MemoryApp.h"
    #include "TimerUtilityVersion.h"
#endif
#ifdef AMP_USE_TRILINOS
    #include "Trilinos_version.h"
#endif
#ifdef AMP_USE_LIBMESH
    #include "libmesh/libmesh_version.h"
#endif
#ifdef AMP_USE_HDF5
    #include "H5public.h"
#endif
#ifdef AMP_USE_SUNDIALS
    #include "sundials/sundials_config.h"
#endif
#ifdef AMP_USE_SILO
    #include "silo.h"
#endif
#ifdef AMP_USE_SAMRAI
    #undef NULL_USE
    #include "SAMRAI/tbox/Logger.h"
    #include "SAMRAI/tbox/SAMRAIManager.h"
    #include "SAMRAI/tbox/StartupShutdownManager.h"
#endif
#ifdef AMP_USE_HYPRE
    #include "HYPRE_config.h"
#endif
#ifdef AMP_USE_KOKKOS
    #include "AMP/utils/KokkosManager.h"
#endif
// clang-format on

#include <algorithm>
#include <array>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <new>
#include <sstream>
#include <stdexcept>
#include <vector>


#define NULL_USE( variable )                    \
    do {                                        \
        if ( 0 ) {                              \
            auto static t = (char *) &variable; \
            t++;                                \
        }                                       \
    } while ( 0 )


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
        Utilities::sleep_ms( ( 100 * comm.getRank() ) / comm.getSize() );
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
static void terminate_AMP2( StackTrace::abort_error &err )
{
    printed_stack = true;
    StackTrace::Utilities::terminate( err );
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
#ifdef AMP_USE_PETSC
    #if PETSC_VERSION_LT( 3, 7, 5 )
        #error AMP only supports PETSc 3.7.5 or greater
    #endif
PetscErrorCode petsc_err_handler( MPI_Comm,
                                  int line,
                                  const char *dir,
                                  const char *file,
                                  PetscErrorCode,
                                  PetscErrorType,
                                  const char *buf,
                                  void * )
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
 * Functions to count resources                                              *
 ****************************************************************************/
static std::map<std::string, std::pair<int, int>> resourceMap;
void AMPManager::incrementResource( const std::string &resource )
{
    auto it = resourceMap.find( resource );
    if ( it == resourceMap.end() )
        resourceMap.insert( std::make_pair( resource, std::pair<int, int>( 1, 0 ) ) );
    else
        it->second.first++;
}
void AMPManager::decrementResource( const std::string &resource )
{
    auto it = resourceMap.find( resource );
    if ( it == resourceMap.end() )
        resourceMap.insert( std::make_pair( resource, std::pair<int, int>( 0, 1 ) ) );
    else
        it->second.second++;
}


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
    StackTrace::Utilities::setAbortBehavior( !AMPManager::use_MPI_Abort, 3 );
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
    // Initialze call stack
    if ( comm_world.getSize() == 1 )
        abort_stackType = std::min( abort_stackType, 2 );
    if ( abort_stackType == 3 )
        StackTrace::globalCallStackInitialize( comm_world.getCommunicator() );
    StackTrace::setDefaultStackType( static_cast<StackTrace::printStackType>( abort_stackType ) );
    // Initialize the random number generator
    AMP::RNG::initialize( 123 );
    // Initialize cuda
    start_CUDA();
    // Initialize Kokkos
    start_Kokkos( argc, argv );
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
        AMP_ERROR( "AMP is not initialized, did you forget to call startup" );
    if ( initialized == -1 )
        AMP_ERROR( "Calling shutdown more than once is invalid" );
    initialized = -1;
    // Clear error handlers
    clearHandlers();
    // Disable MPI_Abort
    AMPManager::use_MPI_Abort = false;
    StackTrace::Utilities::setAbortBehavior( true, 2 );
    // Shutdown the parallel IO
    stopLogging();
    // Shutdown the profiler
    PROFILE_DISABLE();
    // Syncronize all ranks
    comm_world.barrier();
    // shutdown Kokkos
    stop_Kokkos();
    // Shutdown SAMRAI
    double SAMRAI_time = stop_SAMRAI();
    // Shudown PETSc
    double petsc_time = stop_PETSc();
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
    // Print resource counts
    for ( const auto &x : resourceMap ) {
        if ( x.second.first != x.second.second ) {
            printf( "%s:\n", x.first.data() );
            printf( "    %i created\n", x.second.first );
            printf( "    %i destroyed\n", x.second.second );
        }
    }
    resourceMap.clear();
    // Shutdown timer and print memory leaks on rank 0
    PROFILE_DISABLE();
#ifdef AMP_USE_TIMER
    auto memory = MemoryApp::getMemoryStats();
    if ( rank == 0 && memory.N_new > memory.N_delete )
        MemoryApp::print( std::cout );
#endif
    // Wait 50 milli-seconds for all processors to finish
    Utilities::sleep_ms( 50 );
}


/****************************************************************************
 * Restart  the AMP package                                                  *
 ****************************************************************************/
void AMPManager::restart()
{
    if ( initialized != 1 )
        AMP_ERROR( "AMP is not initialized or has been shutdown" );
        // Restart SAMRAI
#ifdef AMP_USE_SAMRAI
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::startup();
#endif
}


/****************************************************************************
 * Function to start/stop SAMRAI                                             *
 ****************************************************************************/
double AMPManager::start_SAMRAI()
{
    double time = 0;
#ifdef AMP_USE_SAMRAI
    double start = Utilities::time();
    #ifdef AMP_USE_MPI
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
#ifdef AMP_USE_SAMRAI
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
#ifdef AMP_USE_SAMRAI
class SAMRAIAbortAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    void logMessage( const std::string &msg, const std::string &file, const int line ) override
    {
        AMP::Utilities::abort( msg, file, line );
    }
    SAMRAIAbortAppender()           = default;
    ~SAMRAIAbortAppender() override = default;
};
#endif


/****************************************************************************
 * Function to start/stop PETSc                                              *
 ****************************************************************************/
double AMPManager::start_PETSc()
{
    double time = 0;
#ifdef AMP_USE_PETSC
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
    #ifndef AMP_USE_MPI
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
#ifdef AMP_USE_PETSC
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
static inline void addArg( const std::string &arg, std::vector<char *> &args )
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
    if ( !properties.initialize_CUDA )
        return 0;
    if ( properties.bind_process_to_accelerator ) {
        auto nodeComm = comm_world.splitByNode();
        auto nodeRank = nodeComm.getRank();
        int deviceCount;
        cudaGetDeviceCount( &deviceCount ); // How many GPUs?
        int device_id = nodeRank % deviceCount;
        cudaSetDevice( device_id ); // Map MPI-process to a GPU
    }

    void *tmp;
    cudaMallocManaged( &tmp, 10, cudaMemAttachGlobal );
    cudaFree( tmp );
#endif
    return 0;
}


/****************************************************************************
 * Function to start/stop Kokkos                                             *
 ****************************************************************************/

double AMPManager::start_Kokkos( int argc, char **argv )
{
#ifdef AMP_USE_KOKKOS
    AMP::Utilities::initializeKokkos( argc, argv );
#else
    NULL_USE( argc );
    NULL_USE( argv );
#endif
    return 0;
}

double AMPManager::stop_Kokkos()
{
#ifdef AMP_USE_KOKKOS
    AMP::Utilities::finalizeKokkos();
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
#ifdef AMP_USE_PETSC
    PetscPopSignalHandler();
    PetscPopErrorHandler();
    PetscPushErrorHandler( &petsc_err_handler, PETSC_NULL );
#endif
// Set the error handlers for SAMRAI
#ifdef AMP_USE_SAMRAI
    SAMRAI::tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit( true );
    SAMRAI::tbox::SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort( true );
    auto appender = std::make_shared<SAMRAIAbortAppender>();
    SAMRAI::tbox::Logger::getInstance()->setAbortAppender( appender );
#endif
    // Set the terminate routine for runtime errors
    StackTrace::Utilities::setErrorHandlers( terminate_AMP2 );
    // Set atexit function
    std::atexit( exitFun );
    int err = std::at_quick_exit( exitFun );
    AMP_ASSERT( err == 0 );
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
#ifdef AMP_USE_PETSC
    PetscPopSignalHandler();
    PetscPopErrorHandler();
#endif
}


/****************************************************************************
 *  Functions to handle MPI errors                                           *
 ****************************************************************************/
void AMPManager::setMPIErrorHandler()
{
#ifdef AMP_USE_MPI
    StackTrace::setMPIErrorHandler( comm_world.getCommunicator() );
    #ifdef AMP_USE_SAMRAI
    StackTrace::setMPIErrorHandler( SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getCommunicator() );
    #endif
#endif
}
void AMPManager::clearMPIErrorHandler()
{
#ifdef AMP_USE_MPI
    StackTrace::clearMPIErrorHandler( comm_world.getCommunicator() );
    #ifdef AMP_USE_SAMRAI
    StackTrace::clearMPIErrorHandler(
        SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getCommunicator() );
    #endif
#endif
}


/****************************************************************************
 * Empty constructor to setup default AMPManagerProperties                   *
 ****************************************************************************/
AMPManagerProperties::AMPManagerProperties() : COMM_WORLD( AMP_COMM_WORLD ) {}


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
#ifdef AMP_USE_TIMER
    out << "ProfilerApp: " << TIMER_VERSION << std::endl;
#endif
#ifdef AMP_USE_SAMRAI
    out << "SAMRAI: " << SAMRAI_VERSION_MAJOR << "." << SAMRAI_VERSION_MINOR << "."
        << SAMRAI_VERSION_PATCHLEVEL << std::endl;
#endif
#ifdef AMP_USE_PETSC
    out << "PETSc: " << PETSC_VERSION_MAJOR << "." << PETSC_VERSION_MINOR << "."
        << PETSC_VERSION_SUBMINOR << std::endl;
#endif
#ifdef AMP_USE_TRILINOS
    out << "Trilinos: " << TRILINOS_VERSION_STRING << std::endl;
#endif
#ifdef AMP_USE_SUNDIALS
    #ifdef SUNDIALS_PACKAGE_VERSION
    out << "Sundials: " << SUNDIALS_PACKAGE_VERSION << std::endl;
    #elif defined( SUNDIALS_VERSION )
    out << "Sundials: " << SUNDIALS_VERSION << std::endl;
    #endif
#endif
#ifdef HYPRE_RELEASE_VERSION
    out << "Hypre: " << HYPRE_RELEASE_VERSION << std::endl;
#elif defined( HYPRE_PACKAGE_VERSION )
    out << "Hypre: " << HYPRE_PACKAGE_VERSION << std::endl;
#endif
#ifdef AMP_USE_LIBMESH
    out << "libMesh: " << libMesh::get_libmesh_version() << std::endl;
#endif
#ifdef AMP_USE_HDF5
    out << "HDF5: " << H5_VERS_MAJOR << "." << H5_VERS_MINOR << "." << H5_VERS_RELEASE << std::endl;
#endif
#ifdef AMP_USE_SILO
    out << "SILO: " << SILO_VERS_MAJ << "." << SILO_VERS_MIN << "." << SILO_VERS_PAT << std::endl;
#endif
#ifdef AMP_USE_MPI
    out << "MPI: " << AMP::AMP_MPI::info();
#endif
#ifdef AMP_USE_LAPACK_WRAPPERS
    out << "Lapack: " << Lapack<double>::info();
#endif
    return out.str();
}


} // namespace AMP
