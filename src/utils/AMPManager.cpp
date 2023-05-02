#include "AMP/utils/AMPManager.h"
#include "AMP/AMP_TPLs.h"
#include "AMP/IO/PIO.h"
#include "AMP/operators/OperatorFactory.h"
#include "AMP/solvers/SolverFactory.h"
#include "AMP/time_integrators/TimeIntegratorFactory.h"
#include "AMP/utils/AMP_MPI.I"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/FactoryStrategy.hpp"
#include "AMP/utils/Utilities.h"

#include "ProfilerApp.h"
#include "StackTrace/ErrorHandlers.h"
#include "StackTrace/StackTrace.h"
#include "StackTrace/Utilities.h"


// Include external packages for startup/shutdown
// clang-format off
#undef NULL_USE
#ifdef USE_CUDA
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include "AMP/utils/cuda/helper_cuda.h"
#endif
#ifdef AMP_USE_PETSC
    #include "petsc.h"
    #include "petscerror.h"
#endif
#ifdef AMP_USE_TIMER
    #include "MemoryApp.h"
#endif
#ifdef AMP_USE_SAMRAI
    #include "SAMRAI/tbox/Logger.h"
    #include "SAMRAI/tbox/SAMRAIManager.h"
    #include "SAMRAI/tbox/Schedule.h"
    #include "SAMRAI/tbox/StartupShutdownManager.h"
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


// Forward declares
namespace AMP {
class KeyData;
}
namespace AMP::Materials {
class Material;
}


namespace AMP {


// Initialize static member variables
int AMPManager::initialized                 = 0;
int AMPManager::abort_stackType             = 3;
bool AMPManager::use_MPI_Abort              = true;
bool AMPManager::print_times                = false;
int AMPManager::argc                        = 0;
char **AMPManager::argv                     = nullptr;
AMPManagerProperties AMPManager::properties = AMPManagerProperties();


/****************************************************************************
 *  Get the global communicator                                              *
 ****************************************************************************/
static AMP_MPI comm_world = AMP::AMP_MPI( MPI_COMM_NULL );
const AMP_MPI &AMPManager::getCommWorld() { return comm_world; }
void AMPManager::setCommWorld( const AMP::AMP_MPI &comm ) { comm_world = comm; }


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
static size_t getStartupMemoryAllocations()
{
#ifdef AMP_USE_TIMER
    auto memory = MemoryApp::getMemoryStats();
    if ( memory.N_new > memory.N_delete )
        return memory.N_new - memory.N_delete;
#endif
    return 0;
}

static size_t N_memory_startup = getStartupMemoryAllocations();
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
    comm_world = AMP_MPI( AMP_COMM_WORLD );
    if ( properties.COMM_WORLD != AMP_COMM_WORLD )
        comm_world = AMP_MPI( properties.COMM_WORLD );
    // Initialize cuda
    start_CUDA();
    // Initialize Kokkos
    start_Kokkos( argc, argv );
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
    // Set the signal/terminate handlers
    StackTrace::Utilities::setErrorHandlers();
    setHandlers();
    // Initialization finished
    initialized = 1;
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
        if ( SAMRAI_time != 0 )
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
    // Shutdown SAMRAI
    double SAMRAI_time = stop_SAMRAI();
    // Shudown PETSc
    double petsc_time = stop_PETSc();
    // shutdown Kokkos
    stop_Kokkos();
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
    // Clear the factories
    AMP::FactoryStrategy<AMP::KeyData>::clear();
    AMP::FactoryStrategy<AMP::Materials::Material>::clear();
    AMP::Operator::OperatorFactory::clear();
    AMP::Solver::SolverFactory::clear();
    AMP::TimeIntegrator::TimeIntegratorFactory::clear();
    // Shutdown timer and print memory leaks on rank 0
    PROFILE_DISABLE();
#ifdef AMP_USE_TIMER
    auto memory = MemoryApp::getMemoryStats();
    if ( rank == 0 && memory.N_new > ( memory.N_delete + N_memory_startup ) ) {
        std::cout << std::endl << "N_memory_startup: " << N_memory_startup << std::endl;
        MemoryApp::print( std::cout );
    }
#endif
    // Wait 50 milli-seconds for all processors to finish
    Utilities::sleep_ms( 50 );
}


/****************************************************************************
 * Restart AMP                                                               *
 ****************************************************************************/
void AMPManager::restart()
{
    if ( initialized != 1 )
        AMP_ERROR( "AMP is not initialized or has been shutdown" );
#ifdef AMP_USE_SAMRAI
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::startup();
#endif
}


/****************************************************************************
 * Function to start/stop SAMRAI                                             *
 ****************************************************************************/
#ifdef AMP_USE_SAMRAI
template<typename T>
class hasClearTimers
{
private:
    template<typename C>
    static char &test( decltype( &C::clearTimers ) );
    template<typename C>
    static int &test( ... );

public:
    static constexpr bool value = sizeof( test<T>( 0 ) ) == sizeof( char );
};
template<typename T>
typename std::enable_if_t<hasClearTimers<T>::value, void> clearTimers( const T &obj )
{
    obj.clearTimers();
}
template<typename T>
typename std::enable_if_t<!hasClearTimers<T>::value, void> clearTimers( const T & )
{
}
double AMPManager::start_SAMRAI()
{
    double start = Utilities::time();
    #ifdef AMP_USE_MPI
    SAMRAI::tbox::SAMRAI_MPI::init( AMP_MPI( AMP_COMM_WORLD ).getCommunicator() );
    #else
    SAMRAI::tbox::SAMRAI_MPI::initMPIDisabled();
    #endif
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();
    SAMRAI::tbox::SAMRAIManager::setMaxNumberPatchDataEntries( 2048 );
    return Utilities::time() - start;
}
double AMPManager::stop_SAMRAI()
{
    double start = Utilities::time();
    SAMRAI::tbox::PIO::finalize();
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();
    clearTimers( SAMRAI::tbox::Schedule() );
    return Utilities::time() - start;
}
#else
double AMPManager::start_SAMRAI() { return 0; }
double AMPManager::stop_SAMRAI() { return 0; }
#endif


/****************************************************************************
 * Function to start/stop PETSc                                              *
 ****************************************************************************/
#ifdef AMP_USE_PETSC
static bool called_PetscInitialize = false;
double AMPManager::start_PETSc()
{
    double time  = 0;
    double start = Utilities::time();
    if ( PetscInitializeCalled ) {
        called_PetscInitialize = false;
    } else {
        int nargsPetsc       = 1;
        const char *noMalloc = "-malloc no";
        char **petscArgs     = const_cast<char **>( &noMalloc );
        PetscInitialize( &nargsPetsc, &petscArgs, PETSC_NULL, PETSC_NULL );
        called_PetscInitialize = true;
    }
    #ifndef AMP_USE_MPI
    // Fix minor bug in petsc where first call to dup returns MPI_COMM_WORLD instead of a new comm
    AMP::AMP_MPI( MPI_COMM_WORLD ).dup();
    #endif
    time = Utilities::time() - start;
    return time;
}
double AMPManager::stop_PETSc()
{
    double time = 0;
    if ( called_PetscInitialize ) {
        double start = Utilities::time();
        PetscPopSignalHandler();
        PetscPopErrorHandler();
        PetscFinalize();
        time = Utilities::time() - start;
    }
    return time;
}
#else
double AMPManager::start_PETSc() { return 0; }
double AMPManager::stop_PETSc() { return 0; }
#endif


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
    checkCudaErrors( cudaMallocManaged( &tmp, 10, cudaMemAttachGlobal ) );
    cudaFree( tmp );
#endif
    return 0;
}


/****************************************************************************
 * Function to start/stop Kokkos                                             *
 ****************************************************************************/
#ifdef AMP_USE_KOKKOS
double AMPManager::start_Kokkos( int argc, char **argv )
{
    AMP::Utilities::initializeKokkos( argc, argv );
    return 0;
}

double AMPManager::stop_Kokkos()
{
    AMP::Utilities::finalizeKokkos();
    return 0;
}
#else
double AMPManager::start_Kokkos( int, char ** ) { return 0; }
double AMPManager::stop_Kokkos() { return 0; }
#endif


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


} // namespace AMP
