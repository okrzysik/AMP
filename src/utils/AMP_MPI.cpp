// This file impliments a wrapper class for MPI functions

// Include AMP headers
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.I"
#include "AMP/utils/Utilities.h"

#include "ProfilerApp.h"
#include "StackTrace/ErrorHandlers.h"
#include "StackTrace/StackTrace.h"

// Include all other headers
#include <algorithm>
#include <chrono>
#include <climits>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <typeinfo>


// Include OS specific headers
#undef USE_WINDOWS
#undef USE_LINUX
#undef USE_MAC
#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )
// We are using windows
#define USE_WINDOWS
#include <process.h>
#include <windows.h>
#define sched_yield() Sleep( 0 )
#elif defined( __APPLE__ )
// Using MAC
#define USE_MAC
#include <sched.h>
#elif defined( __linux ) || defined( __linux__ ) || defined( __unix ) || defined( __posix )
// We are using linux
#define USE_LINUX
#include <sched.h>
#include <unistd.h>
#else
#error Unknown OS
#endif


// Convience defines
#define MPI_CLASS_COMM_NULL AMP_COMM_NULL
#define MPI_CLASS_COMM_SELF AMP_COMM_SELF
#define MPI_CLASS_COMM_WORLD AMP_COMM_WORLD


#if defined( USE_SAMRAI ) && defined( USE_PETSC ) && !defined( USE_MPI )
int MPI_REQUEST_NULL  = 3;
int MPI_ERR_IN_STATUS = 4;
#endif


namespace std {
template<class TYPE>
static inline bool operator<( const std::complex<TYPE> &a, const std::complex<TYPE> &b )
{
    if ( a.real() < b.real() )
        return true;
    if ( a.real() > b.real() )
        return false;
    return a.imag() < b.imag();
}
} // namespace std


namespace AMP {


// Check the alignment
static_assert( sizeof( MPI_CLASS ) % 8 == 0 );


// Define the AMPManager comm_world
AMP_MPI AMPManager::comm_world = AMP::AMP_MPI();


// Initialized the static member variables
volatile uint32_t MPI_CLASS::N_MPI_Comm_created   = 0;
volatile uint32_t MPI_CLASS::N_MPI_Comm_destroyed = 0;
short MPI_CLASS::profile_level                    = 127;


// Default MPI max tag (will be overridden if MPI_TAG_UB attribute is valid)
static const int mpi_max_tag = 0x3FFFFFFF;


// Static data for asyncronous communication without MPI
// Note: these routines may not be thread-safe yet
#ifndef USE_MPI
struct Isendrecv_struct {
    const void *data; // Pointer to data
    int bytes;        // Number of bytes in the message
    int status;       // Status: 1-sending, 2-recieving
    MPI_Comm comm;    // Communicator
    int tag;          // Tag
};
std::map<MPI_Request, Isendrecv_struct> global_isendrecv_list;
static MPI_Request getRequest( MPI_Comm comm, int tag )
{
    MPI_CLASS_ASSERT( tag >= 0 && tag <= mpi_max_tag );
    // Use hashing function: 2^64*0.5*(sqrt(5)-1)
    uint64_t a    = static_cast<uint8_t>( comm ) * 0x9E3779B97F4A7C15;
    uint64_t b    = static_cast<uint8_t>( tag ) * 0x9E3779B97F4A7C15;
    uint64_t hash = a ^ b;
    MPI_Request request;
    memcpy( &request, &hash, sizeof( MPI_Request ) );
    return request;
}
#endif


// Check the mpi error code
#ifdef USE_MPI
inline void check_MPI( int error )
{
    if ( error != MPI_SUCCESS )
        MPI_CLASS_ERROR( "Error calling MPI routine" );
}
#endif


/************************************************************************
 *  Get the MPI version                                                  *
 ************************************************************************/
std::array<int, 2> MPI_CLASS::version()
{
#ifdef USE_MPI
    int MPI_version;
    int MPI_subversion;
    MPI_Get_version( &MPI_version, &MPI_subversion );
    return { MPI_version, MPI_subversion };
#else
    return { 0, 0 };
#endif
}
std::string MPI_CLASS::info()
{
    std::string MPI_info;
#ifdef USE_MPI
#if MPI_VERSION >= 3
    int MPI_version_length = 0;
    char MPI_version_string[MPI_MAX_LIBRARY_VERSION_STRING];
    MPI_Get_library_version( MPI_version_string, &MPI_version_length );
    if ( MPI_version_length > 0 ) {
        MPI_info   = std::string( MPI_version_string, MPI_version_length );
        size_t pos = MPI_info.find( '\n' );
        while ( pos != std::string::npos ) {
            MPI_info.insert( pos + 1, "   " );
            pos = MPI_info.find( '\n', pos + 1 );
        }
    }
#else
    auto tmp = version();
    MPI_info = std::to_string( tmp[0] ) + "." + std::to_string( tmp[0] );
#endif
    size_t pos = MPI_info.find( "\n\n" );
    while ( pos != std::string::npos ) {
        MPI_info.erase( pos + 1, 1 );
        pos = MPI_info.find( '\n', pos + 1 );
    }
    if ( MPI_info.back() == '\n' )
        MPI_info.pop_back();
#endif
    return MPI_info;
}


/************************************************************************
 *  Functions to get/set the process affinities                          *
 ************************************************************************/
int MPI_CLASS::getNumberOfProcessors() { return std::thread::hardware_concurrency(); }
std::vector<int> MPI_CLASS::getProcessAffinity()
{
    std::vector<int> procs;
#ifdef USE_LINUX
    cpu_set_t mask;
    int error = sched_getaffinity( getpid(), sizeof( cpu_set_t ), &mask );
    if ( error != 0 )
        MPI_CLASS_ERROR( "Error getting process affinity" );
    for ( int i = 0; i < (int) sizeof( cpu_set_t ) * CHAR_BIT; i++ ) {
        if ( CPU_ISSET( i, &mask ) )
            procs.push_back( i );
    }
#elif defined( USE_MAC )
    // MAC does not support getting or setting the affinity
    printf( "Warning: MAC does not support getting the process affinity\n" );
    procs.clear();
#elif defined( USE_WINDOWS )
    HANDLE hProc = GetCurrentProcess();
    size_t procMask;
    size_t sysMask;
    PDWORD_PTR procMaskPtr = reinterpret_cast<PDWORD_PTR>( &procMask );
    PDWORD_PTR sysMaskPtr  = reinterpret_cast<PDWORD_PTR>( &sysMask );
    GetProcessAffinityMask( hProc, procMaskPtr, sysMaskPtr );
    for ( int i = 0; i < (int) sizeof( size_t ) * CHAR_BIT; i++ ) {
        if ( ( procMask & 0x1 ) != 0 )
            procs.push_back( i );
        procMask >>= 1;
    }
#else
#error Unknown OS
#endif
    return procs;
}
void MPI_CLASS::setProcessAffinity( const std::vector<int> &procs )
{
#ifdef USE_LINUX
    cpu_set_t mask;
    CPU_ZERO( &mask );
    for ( auto cpu : procs )
        CPU_SET( cpu, &mask );
    int error = sched_setaffinity( getpid(), sizeof( cpu_set_t ), &mask );
    if ( error != 0 )
        MPI_CLASS_ERROR( "Error setting process affinity" );
#elif defined( USE_MAC )
    // MAC does not support getting or setting the affinity
    NULL_USE( procs );
#elif defined( USE_WINDOWS )
    DWORD mask = 0;
    for ( size_t i = 0; i < procs.size(); i++ )
        mask |= ( (DWORD) 1 ) << procs[i];
    HANDLE hProc = GetCurrentProcess();
    SetProcessAffinityMask( hProc, mask );
#else
#error Unknown OS
#endif
}


/************************************************************************
 *  Function to check if MPI is active                                   *
 ************************************************************************/
bool MPI_CLASS::MPI_active()
{
#ifdef USE_MPI
    int initialized = 0, finalized = 0;
    MPI_Initialized( &initialized );
    MPI_Finalized( &finalized );
    return initialized != 0 && finalized == 0;
#else
    return true;
#endif
}
MPI_CLASS::ThreadSupport MPI_CLASS::queryThreadSupport()
{
#ifdef USE_MPI
    int provided = 0;
    MPI_Query_thread( &provided );
    if ( provided == MPI_THREAD_SINGLE )
        return ThreadSupport::SINGLE;
    if ( provided == MPI_THREAD_FUNNELED )
        return ThreadSupport::FUNNELED;
    if ( provided == MPI_THREAD_SERIALIZED )
        return ThreadSupport::SERIALIZED;
    if ( provided == MPI_THREAD_MULTIPLE )
        return ThreadSupport::MULTIPLE;
    return ThreadSupport::SINGLE;
#else
    return ThreadSupport::MULTIPLE;
#endif
}


/************************************************************************
 *  Function to perform a load balance of the given processes            *
 ************************************************************************/
void MPI_CLASS::balanceProcesses( const MPI_CLASS &globalComm,
                                  const int method,
                                  const std::vector<int> &procs,
                                  const int N_min_in,
                                  const int N_max_in )
{
    // Build the list of processors to use
    std::vector<int> cpus = procs;
    if ( cpus.empty() ) {
        for ( int i = 0; i < getNumberOfProcessors(); i++ )
            cpus.push_back( i );
    }
    // Handle the "easy cases"
    if ( method == 1 ) {
        // Trivial case where we do not need any communication
        setProcessAffinity( cpus );
        return;
    }
    // Get the sub-communicator for the current node
    MPI_CLASS nodeComm = globalComm.splitByNode();
    int N_min          = std::min<int>( std::max<int>( N_min_in, 1 ), cpus.size() );
    int N_max          = N_max_in;
    if ( N_max == -1 )
        N_max = cpus.size();
    N_max = std::min<int>( N_max, cpus.size() );
    MPI_CLASS_ASSERT( N_max >= N_min );
    // Perform the load balance within the node
    if ( method == 2 ) {
        int N_proc = cpus.size() / nodeComm.getSize();
        N_proc     = std::max<int>( N_proc, N_min );
        N_proc     = std::min<int>( N_proc, N_max );
        std::vector<int> cpus2( N_proc, -1 );
        for ( int i = 0; i < N_proc; i++ )
            cpus2[i] = cpus[( nodeComm.getRank() * N_proc + i ) % cpus.size()];
        setProcessAffinity( cpus2 );
    } else {
        MPI_CLASS_ERROR( "Unknown method for load balance" );
    }
}


/************************************************************************
 *  Empty constructor                                                    *
 ************************************************************************/
MPI_CLASS::MPI_CLASS()
    : communicator( MPI_COMM_NULL ),
      d_isNull( true ),
      d_manage( false ),
      d_call_abort( true ),
      comm_rank( 0 ),
      comm_size( 1 ),
      d_maxTag( mpi_max_tag ),
      d_currentTag( nullptr ),
      d_ranks( nullptr ),
      d_count( nullptr )
{
}


/************************************************************************
 *  Empty deconstructor                                                  *
 ************************************************************************/
MPI_CLASS::~MPI_CLASS() { reset(); }
void MPI_CLASS::reset()
{
    // Decrement the count if used
    int count = -1;
    if ( d_count != nullptr )
        count = --( *d_count );
    if ( count == 0 ) {
        // We are holding that last reference to the MPI_Comm object, we need to free it
        if ( d_manage ) {
#if defined( USE_MPI ) || defined( USE_PETSC )
            MPI_Comm_set_errhandler( communicator, MPI_ERRORS_ARE_FATAL );
            int err = MPI_Comm_free( &communicator );
            if ( err != MPI_SUCCESS )
                MPI_CLASS_ERROR( "Problem free'ing MPI_Comm object" );
            communicator = MPI_CLASS_COMM_NULL;
            ++N_MPI_Comm_destroyed;
#endif
        }
        if ( d_ranks != nullptr )
            delete[] d_ranks;
        delete d_count;
    }
    if ( d_currentTag == nullptr ) {
        // No tag index
    } else if ( d_currentTag[1] > 1 ) {
        --( d_currentTag[1] );
    } else {
        delete[] d_currentTag;
    }
    d_manage     = false;
    d_count      = nullptr;
    d_ranks      = nullptr;
    comm_rank    = 0;
    comm_size    = 1;
    d_maxTag     = 0;
    d_isNull     = true;
    d_currentTag = nullptr;
    d_call_abort = true;
}


/************************************************************************
 *  Copy constructors                                                    *
 ************************************************************************/
MPI_CLASS::MPI_CLASS( const MPI_CLASS &comm )
    : communicator( comm.communicator ),
      d_isNull( comm.d_isNull ),
      d_manage( comm.d_manage ),
      d_call_abort( comm.d_call_abort ),
      comm_rank( comm.comm_rank ),
      comm_size( comm.comm_size ),
      d_maxTag( comm.d_maxTag ),
      d_currentTag( comm.d_currentTag ),
      d_ranks( comm.d_ranks ),
      d_count( nullptr )
{
    // Initialize the data members to the existing comm object
    if ( d_currentTag != nullptr )
        ++d_currentTag[1];
    d_call_abort = comm.d_call_abort;
    // Set and increment the count
    d_count = comm.d_count;
    if ( d_count != nullptr )
        ++( *d_count );
}
MPI_CLASS::MPI_CLASS( MPI_CLASS &&rhs ) : MPI_CLASS()
{
    std::swap( communicator, rhs.communicator );
    std::swap( d_isNull, rhs.d_isNull );
    std::swap( d_manage, rhs.d_manage );
    std::swap( d_call_abort, rhs.d_call_abort );
    std::swap( profile_level, rhs.profile_level );
    std::swap( comm_rank, rhs.comm_rank );
    std::swap( comm_size, rhs.comm_size );
    std::swap( d_ranks, rhs.d_ranks );
    std::swap( d_maxTag, rhs.d_maxTag );
    std::swap( d_currentTag, rhs.d_currentTag );
    std::swap( d_count, rhs.d_count );
}


/************************************************************************
 *  Assignment operators                                                 *
 ************************************************************************/
MPI_CLASS &MPI_CLASS::operator=( const MPI_CLASS &comm )
{
    if ( this == &comm ) // protect against invalid self-assignment
        return *this;
    // Destroy the previous object
    this->reset();
    // Initialize the data members to the existing object
    this->communicator = comm.communicator;
    this->comm_rank    = comm.comm_rank;
    this->comm_size    = comm.comm_size;
    this->d_ranks      = comm.d_ranks;
    this->d_isNull     = comm.d_isNull;
    this->d_manage     = comm.d_manage;
    this->d_maxTag     = comm.d_maxTag;
    this->d_call_abort = comm.d_call_abort;
    this->d_currentTag = comm.d_currentTag;
    if ( this->d_currentTag != nullptr )
        ++( this->d_currentTag[1] );
    // Set and increment the count
    this->d_count = comm.d_count;
    if ( this->d_count != nullptr )
        ++( *d_count );
    return *this;
}
MPI_CLASS &MPI_CLASS::operator=( MPI_CLASS &&rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    std::swap( communicator, rhs.communicator );
    std::swap( d_isNull, rhs.d_isNull );
    std::swap( d_manage, rhs.d_manage );
    std::swap( d_call_abort, rhs.d_call_abort );
    std::swap( profile_level, rhs.profile_level );
    std::swap( comm_rank, rhs.comm_rank );
    std::swap( comm_size, rhs.comm_size );
    std::swap( d_ranks, rhs.d_ranks );
    std::swap( d_maxTag, rhs.d_maxTag );
    std::swap( d_currentTag, rhs.d_currentTag );
    std::swap( d_count, rhs.d_count );
    return *this;
}


/************************************************************************
 *  Constructor from existing MPI communicator                           *
 ************************************************************************/
int d_global_currentTag_world1[2] = { 1, 1 };
int d_global_currentTag_world2[2] = { 1, 1 };
int d_global_currentTag_self[2]   = { 1, 1 };
#ifdef USE_MPI
std::atomic_int d_global_count_world1 = { 1 };
std::atomic_int d_global_count_world2 = { 1 };
std::atomic_int d_global_count_self   = { 1 };
#endif
MPI_CLASS::MPI_CLASS( MPI_Comm comm, bool manage )
{
    d_count  = nullptr;
    d_ranks  = nullptr;
    d_manage = false;
    // Check if we are using our version of comm_world
    if ( comm == MPI_CLASS_COMM_WORLD ) {
        communicator = AMP::AMPManager::comm_world.communicator;
    } else if ( comm == MPI_CLASS_COMM_SELF ) {
        communicator = MPI_COMM_SELF;
    } else if ( comm == MPI_CLASS_COMM_NULL ) {
        communicator = MPI_COMM_NULL;
    } else {
        communicator = comm;
    }
#ifdef USE_MPI
    // We are using MPI, use the MPI communicator to initialize the data
    if ( communicator != MPI_COMM_NULL ) {
        // Attach the error handler
        StackTrace::setMPIErrorHandler( communicator );
        // Get the communicator properties
        MPI_Comm_rank( communicator, &comm_rank );
        MPI_Comm_size( communicator, &comm_size );
        int flag, *val;
        int ierr = MPI_Comm_get_attr( communicator, MPI_TAG_UB, &val, &flag );
        MPI_CLASS_ASSERT( ierr == MPI_SUCCESS );
        if ( flag == 0 ) {
            d_maxTag = mpi_max_tag; // The tag is not a valid attribute use default value
        } else {
            d_maxTag = *val;
            if ( d_maxTag < 0 ) {
                d_maxTag = 0x7FFFFFFF;
            } // The maximum tag is > a signed int (set to 2^31-1)
            MPI_CLASS_INSIST( d_maxTag >= 0x7FFF, "maximum tag size is < MPI standard" );
        }
    } else {
        comm_rank = 1;
        comm_size = 0;
        d_maxTag  = mpi_max_tag;
    }
    d_isNull = communicator == MPI_COMM_NULL;
    if ( manage && communicator != MPI_COMM_NULL && communicator != MPI_COMM_SELF &&
         communicator != MPI_COMM_WORLD )
        d_manage = true;
    // Create the count (Note: we do not need to worry about thread safety)
    if ( communicator == AMP::AMPManager::comm_world.communicator ) {
        d_count = &d_global_count_world1;
        ++( *d_count );
    } else if ( communicator == MPI_COMM_WORLD ) {
        d_count = &d_global_count_world2;
        ++( *d_count );
    } else if ( communicator == MPI_COMM_SELF ) {
        d_count = &d_global_count_self;
        ++( *d_count );
    } else if ( communicator == MPI_COMM_NULL ) {
        d_count = nullptr;
    } else {
        d_count  = new std::atomic_int;
        *d_count = 1;
    }
    if ( d_manage ) {
        ++N_MPI_Comm_created;
        // StackTrace::multi_stack_info( StackTrace::getCallStack() ).print( std::cout );
    }
    // Create d_ranks
    if ( comm_size > 1 ) {
        d_ranks    = new int[comm_size];
        d_ranks[0] = -1;
    }
#else
    // We are not using MPI, intialize based on the communicator
    NULL_USE( manage );
    comm_rank = 0;
    comm_size = 1;
    d_maxTag  = mpi_max_tag;
    d_isNull  = communicator == MPI_COMM_NULL;
    if ( d_isNull )
        comm_size    = 0;
#endif
    if ( communicator == AMP::AMPManager::comm_world.communicator ) {
        d_currentTag = d_global_currentTag_world1;
        ++( this->d_currentTag[1] );
    } else if ( communicator == MPI_COMM_WORLD ) {
        d_currentTag = d_global_currentTag_world2;
        ++( this->d_currentTag[1] );
    } else if ( communicator == MPI_COMM_SELF ) {
        d_currentTag = d_global_currentTag_self;
        ++( this->d_currentTag[1] );
    } else if ( communicator == MPI_COMM_NULL ) {
        d_currentTag = nullptr;
    } else {
        d_currentTag    = new int[2];
        d_currentTag[0] = ( d_maxTag <= 0x10000 ) ? 1 : 0x1FFF;
        d_currentTag[1] = 1;
    }
    d_call_abort = true;
}


/************************************************************************
 *  Return the ranks of the communicator in the global comm              *
 ************************************************************************/
std::vector<int> MPI_CLASS::globalRanks() const
{
    // Get my global rank if it has not been set
    static int myGlobalRank = -1;
    if ( myGlobalRank == -1 ) {
#ifdef USE_MPI
        if ( MPI_active() )
            MPI_Comm_rank( AMP::AMPManager::comm_world.communicator, &myGlobalRank );
#else
        myGlobalRank = 0;
#endif
    }
    // Check if we are dealing with a serial or null communicator
    if ( comm_size == 1 )
        return std::vector<int>( 1, myGlobalRank );
    if ( d_ranks == nullptr || communicator == MPI_COMM_NULL )
        return std::vector<int>();
    // Fill d_ranks if necessary
    if ( d_ranks[0] == -1 ) {
        if ( communicator == AMP::AMPManager::comm_world.communicator ) {
            for ( int i = 0; i < comm_size; i++ )
                d_ranks[i] = i;
        } else {

            MPI_CLASS_ASSERT( myGlobalRank != -1 );
            this->allGather( myGlobalRank, d_ranks );
        }
    }
    // Return d_ranks
    return std::vector<int>( d_ranks, d_ranks + comm_size );
}


/************************************************************************
 *  Generate a random number                                             *
 ************************************************************************/
size_t MPI_CLASS::rand() const
{
    size_t val = 0;
    if ( getRank() == 0 ) {
        static std::random_device rd;
        static std::mt19937 gen( rd() );
        static std::uniform_int_distribution<size_t> dist;
        val = dist( gen );
    }
    val = bcast( val, 0 );
    return val;
}


/************************************************************************
 *  Intersect two communicators                                          *
 ************************************************************************/
#ifdef USE_MPI
static inline void MPI_Group_free2( MPI_Group *group )
{
    if ( *group != MPI_GROUP_EMPTY ) {
        // MPICH is fine with free'ing an empty group, OpenMPI crashes
        MPI_Group_free( group );
    }
}
MPI_CLASS MPI_CLASS::intersect( const MPI_CLASS &comm1, const MPI_CLASS &comm2 )
{
    MPI_Group group1 = MPI_GROUP_EMPTY, group2 = MPI_GROUP_EMPTY;
    if ( !comm1.isNull() ) {
        MPI_Group_free2( &group1 );
        MPI_Comm_group( comm1.communicator, &group1 );
    }
    if ( !comm2.isNull() ) {
        MPI_Group_free2( &group2 );
        MPI_Comm_group( comm2.communicator, &group2 );
    }
    MPI_Group group12;
    MPI_Group_intersection( group1, group2, &group12 );
    int compare1, compare2;
    MPI_Group_compare( group1, group12, &compare1 );
    MPI_Group_compare( group2, group12, &compare2 );
    MPI_CLASS new_comm( MPI_CLASS_COMM_NULL );
    int size;
    MPI_Group_size( group12, &size );
    if ( compare1 != MPI_UNEQUAL && size != 0 ) {
        // The intersection matches comm1
        new_comm = comm1;
    } else if ( compare2 != MPI_UNEQUAL && size != 0 ) {
        // The intersection matches comm2
        new_comm = comm2;
    } else if ( comm1.isNull() ) {
        // comm1 is null, we can return safely (comm1 is needed for communication)
    } else {
        // The intersection is smaller than comm1 or comm2
        // Check if the new comm is nullptr for all processors
        int max_size = 0;
        MPI_Allreduce( &size, &max_size, 1, MPI_INT, MPI_MAX, comm1.communicator );
        if ( max_size == 0 ) {
            // We are dealing with completely disjoint sets
            new_comm = MPI_CLASS( MPI_CLASS_COMM_NULL, false );
        } else {
            // Create the new comm
            // Note: OpenMPI crashes if the intersection group is EMPTY for any processors
            // We will set it to SELF for the EMPTY processors, then create a nullptr comm later
            if ( group12 == MPI_GROUP_EMPTY ) {
                MPI_Group_free2( &group12 );
                MPI_Comm_group( MPI_COMM_SELF, &group12 );
            }
            MPI_Comm new_MPI_comm;
            MPI_Comm_create( comm1.communicator, group12, &new_MPI_comm );
            if ( size > 0 ) {
                // This is the valid case where we create a new intersection comm
                new_comm = MPI_CLASS( new_MPI_comm, true );
            } else {
                // We actually want a null comm for this communicator
                new_comm = MPI_CLASS( MPI_CLASS_COMM_NULL, false );
                MPI_Comm_free( &new_MPI_comm );
            }
        }
    }
    MPI_Group_free2( &group1 );
    MPI_Group_free2( &group2 );
    MPI_Group_free2( &group12 );
    return new_comm;
}
#else
MPI_CLASS MPI_CLASS::intersect( const MPI_CLASS &comm1, const MPI_CLASS &comm2 )
{
    if ( comm1.isNull() || comm2.isNull() )
        return MPI_CLASS( MPI_CLASS_COMM_NULL, false );
    MPI_CLASS_ASSERT( comm1.comm_size == 1 && comm2.comm_size == 1 );
    return comm1;
}
#endif


/************************************************************************
 *  Split a comm						                                    *
 ************************************************************************/
MPI_CLASS MPI_CLASS::split( int color, int key ) const
{
    if ( d_isNull ) {
        return MPI_CLASS( MPI_CLASS_COMM_NULL );
    } else if ( comm_size == 1 ) {
        if ( color == -1 )
            return MPI_CLASS( MPI_CLASS_COMM_NULL );
        return dup();
    }
    MPI_Comm new_MPI_comm = MPI_CLASS_COMM_NULL;
#ifdef USE_MPI
    // USE MPI to split the communicator
    if ( color == -1 ) {
        check_MPI( MPI_Comm_split( communicator, MPI_UNDEFINED, key, &new_MPI_comm ) );
    } else {
        check_MPI( MPI_Comm_split( communicator, color, key, &new_MPI_comm ) );
    }
#endif
    // Create the new object
    NULL_USE( key );
    MPI_CLASS new_comm( new_MPI_comm, true );
    new_comm.d_call_abort = d_call_abort;
    return new_comm;
}
MPI_CLASS MPI_CLASS::splitByNode( int key ) const
{
    // Check if we are dealing with a single processor (trivial case)
    if ( comm_size == 1 )
        return this->split( 0, 0 );
    // Get the node name
    std::string name = MPI_CLASS::getNodeName();
    // Gather the names from all ranks
    std::vector<std::string> list( comm_size );
    allGather( name, &list[0] );
    // Create the colors
    std::vector<int> color( comm_size, -1 );
    color[0] = 0;
    for ( int i = 1; i < comm_size; i++ ) {
        const std::string tmp1 = list[i];
        for ( int j = 0; j < i; j++ ) {
            const std::string tmp2 = list[j];
            if ( tmp1 == tmp2 ) {
                color[i] = color[j];
                break;
            }
            color[i] = color[i - 1] + 1;
        }
    }
    MPI_CLASS new_comm = this->split( color[comm_rank], key );
    return new_comm;
}


/************************************************************************
 *  Duplicate an exisiting comm object                                   *
 ************************************************************************/
MPI_CLASS MPI_CLASS::dup() const
{
    if ( d_isNull )
        return MPI_CLASS( MPI_CLASS_COMM_NULL );
    MPI_Comm new_MPI_comm = communicator;
#if defined( USE_MPI ) || defined( USE_PETSC )
    // USE MPI to duplicate the communicator
    MPI_Comm_dup( communicator, &new_MPI_comm );
#else
    static MPI_Comm uniqueGlobalComm = 11;
    new_MPI_comm = uniqueGlobalComm;
    uniqueGlobalComm++;
#endif
    // Create the new comm object
    MPI_CLASS new_comm( new_MPI_comm, true );
    new_comm.d_isNull     = d_isNull;
    new_comm.d_call_abort = d_call_abort;
    return new_comm;
}


/************************************************************************
 *  Get the node name                                                    *
 ************************************************************************/
std::string MPI_CLASS::getNodeName()
{
#ifdef USE_MPI
    int length;
    char name[MPI_MAX_PROCESSOR_NAME + 1];
    memset( name, 0, MPI_MAX_PROCESSOR_NAME + 1 );
    MPI_Get_processor_name( name, &length );
    return std::string( name );
#else
    return "Node0";
#endif
}


/************************************************************************
 *  Overload operator ==                                                 *
 ************************************************************************/
bool MPI_CLASS::operator==( const MPI_CLASS &comm ) const
{
    return communicator == comm.communicator;
}


/************************************************************************
 *  Overload operator !=                                                 *
 ************************************************************************/
bool MPI_CLASS::operator!=( const MPI_CLASS &comm ) const
{
    return communicator != comm.communicator;
}


/************************************************************************
 *  Overload operator <                                                  *
 ************************************************************************/
bool MPI_CLASS::operator<( const MPI_CLASS &comm ) const
{
    MPI_CLASS_ASSERT( !this->d_isNull && !comm.d_isNull );
    bool flag = true;
    // First check if either communicator is NULL
    if ( this->d_isNull )
        return false;
    if ( comm.d_isNull )
        flag = false;
    // Use compare to check if the comms are equal
    if ( compare( comm ) != 0 )
        return false;
    // Check that the size of the other communicator is > the current communicator size
    if ( comm_size >= comm.comm_size )
        flag = false;
// Check the union of the communicator groups
// this is < comm iff this group is a subgroup of comm's group
#ifdef USE_MPI
    MPI_Group group1 = MPI_GROUP_EMPTY, group2 = MPI_GROUP_EMPTY, group12 = MPI_GROUP_EMPTY;
    if ( !d_isNull )
        MPI_Comm_group( communicator, &group1 );
    if ( !comm.d_isNull )
        MPI_Comm_group( comm.communicator, &group2 );
    MPI_Group_union( group1, group2, &group12 );
    int compare;
    MPI_Group_compare( group2, group12, &compare );
    if ( compare == MPI_UNEQUAL )
        flag = false;
    MPI_Group_free( &group1 );
    MPI_Group_free( &group2 );
    MPI_Group_free( &group12 );
#endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce( flag );
}


/************************************************************************
 *  Overload operator <=                                                 *
 ************************************************************************/
bool MPI_CLASS::operator<=( const MPI_CLASS &comm ) const
{
    MPI_CLASS_ASSERT( !this->d_isNull && !comm.d_isNull );
    bool flag = true;
    // First check if either communicator is NULL
    if ( this->d_isNull )
        return false;
    if ( comm.d_isNull )
        flag = false;
#ifdef USE_MPI
    int world_size = 0;
    MPI_Comm_size( MPI_COMM_WORLD, &world_size );
    if ( comm.getSize() == world_size )
        return true;
    if ( getSize() == 1 && !comm.d_isNull )
        return true;
#endif
    // Use compare to check if the comms are equal
    if ( compare( comm ) != 0 )
        return true;
    // Check that the size of the other communicator is > the current communicator size
    // this is <= comm iff this group is a subgroup of comm's group
    if ( comm_size > comm.comm_size )
        flag = false;
// Check the unnion of the communicator groups
#ifdef USE_MPI
    MPI_Group group1, group2, group12;
    MPI_Comm_group( communicator, &group1 );
    MPI_Comm_group( comm.communicator, &group2 );
    MPI_Group_union( group1, group2, &group12 );
    int compare;
    MPI_Group_compare( group2, group12, &compare );
    if ( compare == MPI_UNEQUAL )
        flag = false;
    MPI_Group_free( &group1 );
    MPI_Group_free( &group2 );
    MPI_Group_free( &group12 );
#endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce( flag );
}


/************************************************************************
 *  Overload operator >                                                  *
 ************************************************************************/
bool MPI_CLASS::operator>( const MPI_CLASS &comm ) const
{
    bool flag = true;
    // First check if either communicator is NULL
    if ( this->d_isNull )
        return false;
    if ( comm.d_isNull )
        flag = false;
    // Use compare to check if the comms are equal
    if ( compare( comm ) != 0 )
        return false;
    // Check that the size of the other communicator is > the current communicator size
    if ( comm_size <= comm.comm_size )
        flag = false;
// Check the unnion of the communicator groups
// this is > comm iff comm's group is a subgroup of this group
#ifdef USE_MPI
    MPI_Group group1 = MPI_GROUP_EMPTY, group2 = MPI_GROUP_EMPTY, group12 = MPI_GROUP_EMPTY;
    if ( !d_isNull )
        MPI_Comm_group( communicator, &group1 );
    if ( !comm.d_isNull )
        MPI_Comm_group( comm.communicator, &group2 );
    MPI_Group_union( group1, group2, &group12 );
    int compare;
    MPI_Group_compare( group1, group12, &compare );
    if ( compare == MPI_UNEQUAL )
        flag = false;
    MPI_Group_free( &group1 );
    MPI_Group_free( &group2 );
    MPI_Group_free( &group12 );
#else
    NULL_USE( comm );
#endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce( flag );
}


/************************************************************************
 *  Overload operator >=                                                 *
 ************************************************************************/
bool MPI_CLASS::operator>=( const MPI_CLASS &comm ) const
{
    bool flag = true;
    // First check if either communicator is NULL
    if ( this->d_isNull )
        return false;
    if ( comm.d_isNull )
        flag = false;
#ifdef USE_MPI
    int world_size = 0;
    MPI_Comm_size( MPI_COMM_WORLD, &world_size );
    if ( getSize() == world_size )
        return true;
    if ( comm.getSize() == 1 && !comm.d_isNull )
        return true;
#endif
    // Use compare to check if the comms are equal
    if ( compare( comm ) != 0 )
        return true;
    // Check that the size of the other communicator is > the current communicator size
    if ( comm_size < comm.comm_size )
        flag = false;
// Check the unnion of the communicator groups
// this is >= comm iff comm's group is a subgroup of this group
#ifdef USE_MPI
    MPI_Group group1 = MPI_GROUP_EMPTY, group2 = MPI_GROUP_EMPTY, group12 = MPI_GROUP_EMPTY;
    if ( !d_isNull )
        MPI_Comm_group( communicator, &group1 );
    if ( !comm.d_isNull )
        MPI_Comm_group( comm.communicator, &group2 );
    MPI_Group_union( group1, group2, &group12 );
    int compare;
    MPI_Group_compare( group1, group12, &compare );
    if ( compare == MPI_UNEQUAL )
        flag = false;
    MPI_Group_free( &group1 );
    MPI_Group_free( &group2 );
    MPI_Group_free( &group12 );
#endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce( flag );
}


/************************************************************************
 *  Compare two comm objects                                             *
 ************************************************************************/
int MPI_CLASS::compare( const MPI_CLASS &comm ) const
{
    if ( communicator == comm.communicator )
        return 1;
#ifdef USE_MPI
    if ( d_isNull || comm.d_isNull )
        return 0;
    int result;
    check_MPI( MPI_Comm_compare( communicator, comm.communicator, &result ) );
    if ( result == MPI_IDENT )
        return 2;
    else if ( result == MPI_CONGRUENT )
        return 3;
    else if ( result == MPI_SIMILAR )
        return 4;
    else if ( result == MPI_UNEQUAL )
        return 0;
    MPI_CLASS_ERROR( "Unknown results from comm compare" );
#else
    if ( comm.communicator == MPI_COMM_NULL || communicator == MPI_COMM_NULL )
        return 0;
    else
        return 3;
#endif
    return 0;
}


/************************************************************************
 *  Abort the program.                                                   *
 ************************************************************************/
void MPI_CLASS::setCallAbortInSerialInsteadOfExit( bool flag ) { d_call_abort = flag; }
void MPI_CLASS::abort() const
{
#ifdef USE_MPI
    MPI_Comm comm = communicator;
    if ( comm == MPI_COMM_NULL )
        comm = MPI_COMM_WORLD;
    if ( !MPI_active() ) {
        // MPI is not availible
        exit( -1 );
    } else if ( comm_size > 1 ) {
        MPI_Abort( comm, -1 );
    } else if ( d_call_abort ) {
        MPI_Abort( comm, -1 );
    } else {
        exit( -1 );
    }
#else
    exit( -1 );
#endif
}


/************************************************************************
 *  newTag                                                               *
 ************************************************************************/
int MPI_CLASS::newTag()
{
#ifdef USE_MPI
    // Syncronize the processes to ensure all ranks enter this call
    // Needed so the count will match
    barrier();
    // Return and increment the tag
    int tag = ( *d_currentTag )++;
    MPI_CLASS_INSIST( tag <= d_maxTag, "Maximum number of tags exceeded\n" );
    return tag;
#else
    static int globalCurrentTag = 1;
    return globalCurrentTag++;
#endif
}


/************************************************************************
 *  allReduce                                                            *
 ************************************************************************/
static inline void setBit( uint64_t *x, size_t index )
{
    size_t i      = index >> 6;
    size_t j      = index & 0x3F;
    uint64_t mask = ( (uint64_t) 0x01 ) << j;
    x[i]          = x[i] | mask;
}
static inline bool readBit( const uint64_t *x, size_t index )
{
    size_t i      = index >> 6;
    size_t j      = index & 0x3F;
    uint64_t mask = ( (uint64_t) 0x01 ) << j;
    return ( x[i] & mask ) != 0;
}
bool MPI_CLASS::allReduce( const bool value ) const
{
    bool ret = value;
    if ( comm_size > 1 ) {
#ifdef USE_MPI
        MPI_Allreduce(
            (void *) &value, (void *) &ret, 1, MPI_UNSIGNED_CHAR, MPI_MIN, communicator );
#else
        NULL_USE( value );
        MPI_CLASS_ERROR( "This shouldn't be possible" );
#endif
    }
    return ret;
}
void MPI_CLASS::allReduce( std::vector<bool> &x ) const
{
    if ( comm_size <= 1 )
        return;
#ifdef USE_MPI
    size_t N  = ( x.size() + 63 ) / 64;
    auto send = new uint64_t[N];
    auto recv = new uint64_t[N];
    memset( send, 0, N * sizeof( int64_t ) );
    memset( recv, 0, N * sizeof( int64_t ) );
    for ( size_t i = 0; i < x.size(); i++ ) {
        if ( x[i] )
            setBit( send, i );
    }
    MPI_Allreduce( send, recv, N, MPI_UINT64_T, MPI_BAND, communicator );
    for ( size_t i = 0; i < x.size(); i++ )
        x[i] = readBit( recv, i );
    delete[] send;
    delete[] recv;
#else
    NULL_USE( x );
    MPI_CLASS_ERROR( "This shouldn't be possible" );
#endif
}


/************************************************************************
 *  anyReduce                                                            *
 ************************************************************************/
bool MPI_CLASS::anyReduce( const bool value ) const
{
    bool ret = value;
    if ( comm_size > 1 ) {
#ifdef USE_MPI
        MPI_Allreduce(
            (void *) &value, (void *) &ret, 1, MPI_UNSIGNED_CHAR, MPI_MAX, communicator );
#else
        NULL_USE( value );
        MPI_CLASS_ERROR( "This shouldn't be possible" );
#endif
    }
    return ret;
}
void MPI_CLASS::anyReduce( std::vector<bool> &x ) const
{
    if ( comm_size <= 1 )
        return;
#ifdef USE_MPI
    size_t N  = ( x.size() + 63 ) / 64;
    auto send = new uint64_t[N];
    auto recv = new uint64_t[N];
    memset( send, 0, N * sizeof( int64_t ) );
    memset( recv, 0, N * sizeof( int64_t ) );
    for ( size_t i = 0; i < x.size(); i++ ) {
        if ( x[i] )
            setBit( send, i );
    }
    MPI_Allreduce( send, recv, N, MPI_UINT64_T, MPI_BOR, communicator );
    for ( size_t i = 0; i < x.size(); i++ )
        x[i] = readBit( recv, i );
    delete[] send;
    delete[] recv;
#else
    NULL_USE( x );
    MPI_CLASS_ERROR( "This shouldn't be possible" );
#endif
}


/************************************************************************
 *  Perform a global barrier across all processors.                      *
 ************************************************************************/
void MPI_CLASS::barrier() const
{
#ifdef USE_MPI
    MPI_Barrier( communicator );
#endif
}


/************************************************************************
 *  Send/Recv data                                                       *
 *  We need a concrete instantiation of send for use without MPI         *
 ************************************************************************/
#if defined( USE_MPI ) || defined( USE_EXT_MPI )
void MPI_CLASS::sendBytes( const void *buf, int bytes, int recv_proc, int tag ) const
{
    send<char>( (const char *) buf, bytes, recv_proc, tag );
}
void MPI_CLASS::recvBytes( void *buf, int bytes, const int send_proc, int tag ) const
{
    int bytes2 = bytes;
    recv<char>( (char *) buf, bytes2, send_proc, false, tag );
}
MPI_Request MPI_CLASS::IsendBytes( const void *buf, int bytes, int recv_proc, int tag ) const
{
    return Isend<char>( (const char *) buf, bytes, recv_proc, tag );
}
MPI_Request MPI_CLASS::IrecvBytes( void *buf, int bytes, int send_proc, const int tag ) const
{
    return Irecv<char>( (char *) buf, bytes, send_proc, tag );
}
#else
void MPI_CLASS::sendBytes( const void *buf, int bytes, int, int tag ) const
{
    MPI_CLASS_INSIST( tag <= d_maxTag, "Maximum tag value exceeded" );
    MPI_CLASS_INSIST( tag >= 0, "tag must be >= 0" );
    PROFILE_START( "sendBytes", profile_level );
    auto id = getRequest( communicator, tag );
    auto it = global_isendrecv_list.find( id );
    MPI_CLASS_INSIST( it == global_isendrecv_list.end(),
                "send must be paired with a previous call to irecv in serial" );
    MPI_CLASS_ASSERT( it->second.status == 2 );
    memcpy( (char *) it->second.data, buf, bytes );
    global_isendrecv_list.erase( it );
    PROFILE_STOP( "sendBytes", profile_level );
}
void MPI_CLASS::recvBytes( void *buf, int bytes, const int, int tag ) const
{
    MPI_CLASS_INSIST( tag <= d_maxTag, "Maximum tag value exceeded" );
    MPI_CLASS_INSIST( tag >= 0, "tag must be >= 0" );
    PROFILE_START( "recv<char>", profile_level );
    auto id = getRequest( communicator, tag );
    auto it = global_isendrecv_list.find( id );
    MPI_CLASS_INSIST( it != global_isendrecv_list.end(),
                "recv must be paired with a previous call to isend in serial" );
    MPI_CLASS_ASSERT( it->second.status == 1 );
    MPI_CLASS_ASSERT( it->second.bytes == bytes );
    memcpy( buf, it->second.data, bytes );
    global_isendrecv_list.erase( it );
    PROFILE_STOP( "recv<char>", profile_level );
}
MPI_Request MPI_CLASS::IsendBytes( const void *buf, int bytes, int, int tag ) const
{
    MPI_CLASS_INSIST( tag <= d_maxTag, "Maximum tag value exceeded" );
    MPI_CLASS_INSIST( tag >= 0, "tag must be >= 0" );
    PROFILE_START( "IsendBytes", profile_level );
    auto id = getRequest( communicator, tag );
    auto it = global_isendrecv_list.find( id );
    if ( it == global_isendrecv_list.end() ) {
        // We are calling isend first
        Isendrecv_struct data;
        data.bytes  = bytes;
        data.data   = buf;
        data.status = 1;
        data.comm   = communicator;
        data.tag    = tag;
        global_isendrecv_list.insert( std::pair<MPI_Request, Isendrecv_struct>( id, data ) );
    } else {
        // We called irecv first
        MPI_CLASS_ASSERT( it->second.status == 2 );
        MPI_CLASS_ASSERT( it->second.bytes == bytes );
        memcpy( (char *) it->second.data, buf, bytes );
        global_isendrecv_list.erase( it );
    }
    PROFILE_STOP( "IsendBytes", profile_level );
    return id;
}
MPI_Request MPI_CLASS::IrecvBytes( void *buf, const int bytes, const int, const int tag ) const
{
    MPI_CLASS_INSIST( tag <= d_maxTag, "Maximum tag value exceeded" );
    MPI_CLASS_INSIST( tag >= 0, "tag must be >= 0" );
    PROFILE_START( "Irecv<char>", profile_level );
    auto id = getRequest( communicator, tag );
    auto it = global_isendrecv_list.find( id );
    if ( it == global_isendrecv_list.end() ) {
        // We are calling Irecv first
        Isendrecv_struct data;
        data.bytes  = bytes;
        data.data   = buf;
        data.status = 2;
        data.comm   = communicator;
        data.tag    = tag;
        global_isendrecv_list.insert( std::pair<MPI_Request, Isendrecv_struct>( id, data ) );
    } else {
        // We called Isend first
        MPI_CLASS_ASSERT( it->second.status == 1 );
        MPI_CLASS_ASSERT( it->second.bytes == bytes );
        memcpy( buf, it->second.data, bytes );
        global_isendrecv_list.erase( it );
    }
    PROFILE_STOP( "Irecv<char>", profile_level );
    return id;
}
#endif


/************************************************************************
 *  call_allGather                                                       *
 *  Note: these specializations are only called when using MPI.          *
 ************************************************************************/
#ifdef USE_MPI
template<>
void MPI_CLASS::call_allGather<std::string>( const std::string &x_in, std::string *x_out ) const
{
    // Get the bytes recvied per processor
    std::vector<int> recv_cnt( comm_size, 0 );
    allGather<int>( (int) x_in.size() + 1, &recv_cnt[0] );
    std::vector<int> recv_disp( comm_size, 0 );
    for ( int i = 1; i < comm_size; i++ )
        recv_disp[i] = recv_disp[i - 1] + recv_cnt[i - 1];
    // Call the vector form of allGather for the char arrays
    char *recv_data = new char[recv_disp[comm_size - 1] + recv_cnt[comm_size - 1]];
    allGather<char>(
        x_in.c_str(), (int) x_in.size() + 1, recv_data, &recv_cnt[0], &recv_disp[0], true );
    for ( int i = 0; i < comm_size; i++ )
        x_out[i] = std::string( &recv_data[recv_disp[i]] );
    delete[] recv_data;
}
template<>
void MPI_CLASS::call_allGather<std::string>(
    const std::string *, int, std::string *, int *, int * ) const
{
    throw std::logic_error( "Implimentation of allgatherv for std::string not implimented yet" );
}
#endif


/************************************************************************
 *  Communicate ranks for communication                                  *
 ************************************************************************/
std::vector<int> MPI_CLASS::commRanks( const std::vector<int> &ranks ) const
{
#ifdef USE_MPI
    // Get a byte array with the ranks to communicate
    auto data1 = new char[comm_size];
    auto data2 = new char[comm_size];
    memset( data1, 0, comm_size );
    memset( data2, 0, comm_size );
    for ( auto &rank : ranks )
        data1[rank] = 1;
    MPI_Alltoall( data1, 1, MPI_CHAR, data2, 1, MPI_CHAR, communicator );
    int N = 0;
    for ( int i = 0; i < comm_size; i++ )
        N += data2[i];
    std::vector<int> ranks_out;
    ranks_out.reserve( N );
    for ( int i = 0; i < comm_size; i++ ) {
        if ( data2[i] )
            ranks_out.push_back( i );
    }
    delete[] data1;
    delete[] data2;
    return ranks_out;
#else
    return ranks;
#endif
}


/************************************************************************
 *  Wait functions                                                       *
 ************************************************************************/
#ifdef USE_MPI
void MPI_CLASS::wait( MPI_Request request )
{
    PROFILE_START( "wait", profile_level );
    MPI_Status status;
    int flag = 0;
    int err  = MPI_Test( &request, &flag, &status );
    MPI_CLASS_ASSERT( err == MPI_SUCCESS ); // Check that the first call is valid
    while ( !flag ) {
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
        // Check if the request has finished
        MPI_Test( &request, &flag, &status );
    }
    PROFILE_STOP( "wait", profile_level );
}
int MPI_CLASS::waitAny( int count, MPI_Request *request )
{
    if ( count == 0 )
        return -1;
    PROFILE_START( "waitAny", profile_level );
    int index   = -1;
    int flag    = 0;
    auto status = new MPI_Status[count];
    int err     = MPI_Testany( count, request, &index, &flag, status );
    MPI_CLASS_ASSERT( err == MPI_SUCCESS ); // Check that the first call is valid
    while ( !flag ) {
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
        // Check if the request has finished
        MPI_Testany( count, request, &index, &flag, status );
    }
    MPI_CLASS_ASSERT( index >= 0 ); // Check that the index is valid
    delete[] status;
    PROFILE_STOP( "waitAny", profile_level );
    return index;
}
void MPI_CLASS::waitAll( int count, MPI_Request *request )
{
    if ( count == 0 )
        return;
    PROFILE_START( "waitAll", profile_level );
    int flag    = 0;
    auto status = new MPI_Status[count];
    int err     = MPI_Testall( count, request, &flag, status );
    MPI_CLASS_ASSERT( err == MPI_SUCCESS ); // Check that the first call is valid
    while ( !flag ) {
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
        // Check if the request has finished
        MPI_Testall( count, request, &flag, status );
    }
    PROFILE_STOP( "waitAll", profile_level );
    delete[] status;
}
std::vector<int> MPI_CLASS::waitSome( int count, MPI_Request *request )
{
    if ( count == 0 )
        return std::vector<int>();
    PROFILE_START( "waitSome", profile_level );
    std::vector<int> indicies( count, -1 );
    auto *status = new MPI_Status[count];
    int outcount = 0;
    int err      = MPI_Testsome( count, request, &outcount, &indicies[0], status );
    MPI_CLASS_ASSERT( err == MPI_SUCCESS );        // Check that the first call is valid
    MPI_CLASS_ASSERT( outcount != MPI_UNDEFINED ); // Check that the first call is valid
    while ( outcount == 0 ) {
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
        // Check if the request has finished
        MPI_Testsome( count, request, &outcount, &indicies[0], status );
    }
    indicies.resize( outcount );
    delete[] status;
    PROFILE_STOP( "waitSome", profile_level );
    return indicies;
}
#else
void MPI_CLASS::wait( MPI_Request request )
{
    PROFILE_START( "wait", profile_level );
    while ( 1 ) {
        // Check if the request is in our list
        if ( global_isendrecv_list.find( request ) == global_isendrecv_list.end() )
            break;
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
    }
    PROFILE_STOP( "wait", profile_level );
}
int MPI_CLASS::waitAny( int count, MPI_Request *request )
{
    if ( count == 0 )
        return -1;
    PROFILE_START( "waitAny", profile_level );
    int index = 0;
    while ( 1 ) {
        // Check if the request is in our list
        bool found_any = false;
        for ( int i = 0; i < count; i++ ) {
            if ( global_isendrecv_list.find( request[i] ) == global_isendrecv_list.end() ) {
                found_any = true;
                index     = i;
            }
        }
        if ( found_any )
            break;
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
    }
    PROFILE_STOP( "waitAny", profile_level );
    return index;
}
void MPI_CLASS::waitAll( int count, MPI_Request *request )
{
    if ( count == 0 )
        return;
    PROFILE_START( "waitAll", profile_level );
    while ( 1 ) {
        // Check if the request is in our list
        bool found_all = true;
        for ( int i = 0; i < count; i++ ) {
            if ( global_isendrecv_list.find( request[i] ) != global_isendrecv_list.end() )
                found_all = false;
        }
        if ( found_all )
            break;
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
    }
    PROFILE_STOP( "waitAll", profile_level );
}
std::vector<int> MPI_CLASS::waitSome( int count, MPI_Request *request )
{
    if ( count == 0 )
        return std::vector<int>();
    PROFILE_START( "waitSome", profile_level );
    std::vector<int> indicies;
    while ( 1 ) {
        // Check if the request is in our list
        for ( int i = 0; i < count; i++ ) {
            if ( global_isendrecv_list.find( request[i] ) == global_isendrecv_list.end() )
                indicies.push_back( i );
        }
        if ( !indicies.empty() )
            break;
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
    }
    PROFILE_STOP( "waitSome", profile_level );
    return indicies;
}
#endif


/************************************************************************
 *  Probe functions                                                      *
 ************************************************************************/
#ifdef USE_MPI
int MPI_CLASS::Iprobe( int source, int tag ) const
{
    MPI_CLASS_INSIST( tag <= d_maxTag, "Maximum tag value exceeded" );
    MPI_CLASS_INSIST( tag >= 0, "tag must be >= 0" );
    MPI_Status status;
    int flag = 0;
    MPI_Iprobe( source, tag, communicator, &flag, &status );
    if ( flag == 0 )
        return -1;
    int count;
    MPI_Get_count( &status, MPI_BYTE, &count );
    MPI_CLASS_ASSERT( count >= 0 );
    return count;
}
int MPI_CLASS::probe( int source, int tag ) const
{
    MPI_CLASS_INSIST( tag <= d_maxTag, "Maximum tag value exceeded" );
    MPI_CLASS_INSIST( tag >= 0, "tag must be >= 0" );
    MPI_Status status;
    MPI_Probe( source, tag, communicator, &status );
    int count;
    MPI_Get_count( &status, MPI_BYTE, &count );
    MPI_CLASS_ASSERT( count >= 0 );
    return count;
}
#else
int MPI_CLASS::Iprobe( int source, int tag ) const
{
    AMP_ASSERT( source == 0 || source < -1 );
    for ( const auto &tmp : global_isendrecv_list ) {
        const auto &data = tmp.second;
        if ( data.comm == communicator && ( data.tag == tag || tag == -1 ) && data.status == 1 )
            return data.bytes;
    }
    return -1;
}
int MPI_CLASS::probe( int source, int tag ) const
{
    int bytes = Iprobe( source, tag );
    AMP_INSIST( bytes >= 0, "probe called before message started in serial" );
    return bytes;
}
#endif


/************************************************************************
 *  Timer functions                                                      *
 ************************************************************************/
#ifdef USE_MPI
double MPI_CLASS::time() { return MPI_Wtime(); }
double MPI_CLASS::tick() { return MPI_Wtick(); }
#else
double MPI_CLASS::time()
{
    auto t  = std::chrono::system_clock::now();
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>( t.time_since_epoch() );
    return 1e-9 * ns.count();
}
double MPI_CLASS::tick()
{
    auto period = std::chrono::system_clock::period();
    return static_cast<double>( period.num ) / static_cast<double>( period.den );
}
#endif


/************************************************************************
 *  Serialize a block of code across MPI processes                       *
 ************************************************************************/
void MPI_CLASS::serializeStart()
{
#ifdef USE_MPI
    using namespace std::chrono_literals;
    if ( comm_rank == 0 ) {
        // Start rank 0 immediately
    } else {
        // Wait for a message from the previous rank
        MPI_Request request;
        MPI_Status status;
        int flag = false, buf = 0;
        MPI_Irecv( &buf, 1, MPI_INT, comm_rank - 1, 5627, MPI_COMM_WORLD, &request );
        while ( !flag ) {
            MPI_Test( &request, &flag, &status );
            std::this_thread::sleep_for( 50ms );
        }
    }
#endif
}
void MPI_CLASS::serializeStop()
{
#ifdef USE_MPI
    using namespace std::chrono_literals;
    if ( comm_rank < comm_size - 1 ) {
        // Send flag to next rank
        MPI_Send( &comm_rank, 1, MPI_INT, comm_rank + 1, 5627, MPI_COMM_WORLD );
        // Wait for final finished flag
        int flag = false, buf = 0;
        MPI_Request request;
        MPI_Status status;
        MPI_Irecv( &buf, 1, MPI_INT, comm_size - 1, 5627, MPI_COMM_WORLD, &request );
        while ( !flag ) {
            MPI_Test( &request, &flag, &status );
            std::this_thread::sleep_for( 50ms );
        }
    } else {
        // Send final flag to all ranks
        for ( int i = 0; i < comm_size - 1; i++ )
            MPI_Send( &comm_rank, 1, MPI_INT, i, 5627, MPI_COMM_WORLD );
    }
#endif
}


/****************************************************************************
 * Function to start/stop MPI                                                *
 ****************************************************************************/
#ifdef USE_EXT_MPI
static bool called_MPI_Init = false;
#endif
bool MPI_CLASS::MPI_Active()
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
void MPI_CLASS::start_MPI( int argc, char *argv[], int profile_level )
{
    changeProfileLevel( profile_level );
    NULL_USE( argc );
    NULL_USE( argv );
#ifdef USE_EXT_MPI
    if ( MPI_Active() ) {
        called_MPI_Init = false;
    } else {
        int provided;
        int result = MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided );
        if ( result != MPI_SUCCESS )
            MPI_CLASS_ERROR( "AMP was unable to initialize MPI" );
        if ( provided < MPI_THREAD_MULTIPLE )
            AMP::perr << "Warning: Failed to start MPI with MPI_THREAD_MULTIPLE\n";
        called_MPI_Init        = true;
        AMPManager::comm_world = AMP_MPI( MPI_COMM_WORLD );
    }
#else
    AMPManager::comm_world = AMP_MPI( MPI_COMM_SELF );
#endif
}
void MPI_CLASS::stop_MPI()
{
    AMPManager::comm_world = AMP_MPI( AMP_COMM_NULL );
#ifdef USE_EXT_MPI
    int finalized;
    MPI_Finalized( &finalized );
    if ( called_MPI_Init && !finalized ) {
        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Finalize();
        called_MPI_Init = true;
    }
#endif
}


/****************************************************************************
 * call_bcast                                                                *
 ****************************************************************************/
#ifdef USE_EXT_MPI
template<>
void MPI_CLASS::call_bcast<std::string>( std::string *str, int n, int root ) const
{
    // Send the length of the strings
    std::vector<int> length( n, 0 );
    if ( root == comm_rank ) {
        for ( int i = 0; i < n; i++ )
            length[i] = str[i].size();
    }
    bcast( length.data(), n, root );
    // Allocate space for the temporary buffer
    size_t N = 0;
    for ( int i = 0; i < n; i++ )
        N += length[i];
    auto buffer = new char[N];
    // Create and send the buffer
    if ( root == comm_rank ) {
        for ( int i = 0, j = 0; i < n; i++ ) {
            memcpy( &buffer[j], str[i].data(), length[i] );
            j += length[i];
        }
    }
    MPI_Bcast( buffer, N, MPI_CHAR, root, communicator );
    // Unpack the strings
    if ( root != comm_rank ) {
        for ( int i = 0, j = 0; i < n; i++ ) {
            str[i].resize( length[i] );
            memcpy( str[i].data(), &buffer[j], length[i] );
            j += length[i];
        }
    }
    delete[] buffer;
}
#endif


/****************************************************************************
 * Explicit instantiation                                                    *
 ****************************************************************************/
// clang-format off
#define INSTANTIATE( TYPE )                                                             \
    template TYPE MPI_CLASS::sumReduce<TYPE>( const TYPE ) const;                       \
    template void MPI_CLASS::sumReduce<TYPE>( TYPE*, int ) const;                       \
    template void MPI_CLASS::sumReduce<TYPE>( const TYPE*, TYPE*, int ) const;          \
    template TYPE MPI_CLASS::minReduce<TYPE>( const TYPE ) const;                       \
    template void MPI_CLASS::minReduce<TYPE>( TYPE*, int ) const;                       \
    template void MPI_CLASS::minReduce<TYPE>( TYPE*, int, int* ) const;                 \
    template void MPI_CLASS::minReduce<TYPE>( const TYPE*, TYPE*, int, int* ) const;    \
    template TYPE MPI_CLASS::maxReduce<TYPE>( const TYPE ) const;                       \
    template void MPI_CLASS::maxReduce<TYPE>( TYPE*, int ) const;                       \
    template void MPI_CLASS::maxReduce<TYPE>( TYPE*, int, int* ) const;                 \
    template void MPI_CLASS::maxReduce<TYPE>( const TYPE*, TYPE*, int, int* ) const;    \
    template void MPI_CLASS::sumScan<TYPE>( const TYPE*, TYPE*, int ) const;            \
    template void MPI_CLASS::minScan<TYPE>( const TYPE*, TYPE*, int ) const;            \
    template void MPI_CLASS::maxScan<TYPE>( const TYPE*, TYPE*, int ) const;            \
    template TYPE MPI_CLASS::bcast<TYPE>( const TYPE&, int ) const;                     \
    template void MPI_CLASS::bcast<TYPE>( TYPE*, int, int ) const;                      \
    template void MPI_CLASS::send<TYPE>( const TYPE*, int, int, int ) const;            \
    template MPI_Request MPI_CLASS::Isend<TYPE>( const TYPE*, int, int, int ) const;    \
    template void MPI_CLASS::recv<TYPE>( TYPE*, int&, int, const bool, int ) const;     \
    template MPI_Request MPI_CLASS::Irecv<TYPE>( TYPE* buf, int, int, int ) const;      \
    template std::vector<TYPE> MPI_CLASS::allGather<TYPE>( const TYPE & ) const;        \
    template std::vector<TYPE> MPI_CLASS::allGather<TYPE>( const std::vector<TYPE>& ) const; \
    template void MPI_CLASS::allGather<TYPE>( const TYPE&, TYPE* ) const;               \
    template int MPI_CLASS::allGather<TYPE>( const TYPE*, int, TYPE*, int*, int*, bool ) const; \
    template void MPI_CLASS::setGather<TYPE>( std::set<TYPE>& ) const;                  \
    template void MPI_CLASS::allToAll<TYPE>( int, const TYPE*, TYPE* ) const;           \
    template int MPI_CLASS::allToAll<TYPE>( const TYPE*, const int[], const int[], TYPE*, int*, int*, bool ) const;
// Instantiate basic types
INSTANTIATE( char )
INSTANTIATE( int8_t )
INSTANTIATE( uint8_t )
INSTANTIATE( int16_t )
INSTANTIATE( uint16_t )
INSTANTIATE( int32_t )
INSTANTIATE( uint32_t )
INSTANTIATE( int64_t )
INSTANTIATE( uint64_t )
INSTANTIATE( float )
INSTANTIATE( double )
INSTANTIATE( std::complex<float> )
INSTANTIATE( std::complex<double> )
// Instantiate std::string
template std::string MPI_CLASS::bcast<std::string>( const std::string&, int ) const;
template void MPI_CLASS::bcast<std::string>( std::string*, int, int ) const;
template std::vector<std::string> MPI_CLASS::allGather<std::string>( const std::string & ) const;
template void MPI_CLASS::allGather<std::string>( const std::string&, std::string* ) const;
template int MPI_CLASS::allGather<std::string>( const std::string*, int, std::string*, int*, int*, bool ) const;
// clang-format on


} // namespace AMP
