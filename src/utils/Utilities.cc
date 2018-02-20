#define NOMINMAX
#include "AMP/utils/Utilities.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Logger.h"
#include "AMP/utils/PIO.h"

#include "StackTrace/StackTrace.h"

#include <algorithm>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>

#ifdef USE_TIMER
#include "MemoryApp.h"
#endif

// Detect the OS
// clang-format off
#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 ) || defined( _MSC_VER )
    #define USE_WINDOWS
#elif defined( __APPLE__ )
    #define USE_MAC
#elif defined( __linux ) || defined( __unix ) || defined( __posix )
    #define USE_LINUX
    #define USE_NM
#else
    #error Unknown OS
#endif
// clang-format on


// Include system dependent headers
// clang-format off
#ifdef USE_WINDOWS
    #include <process.h>
    #include <psapi.h>
    #include <stdio.h>
    #include <tchar.h>
    #include <windows.h>
#else
    #include <dlfcn.h>
    #include <execinfo.h>
    #include <sched.h>
    #include <sys/time.h>
    #include <ctime>
    #include <unistd.h>
#endif
#ifdef USE_LINUX
    #include <malloc.h>
#endif
#ifdef USE_MAC
    #include <mach/mach.h>
    #include <sys/sysctl.h>
    #include <sys/types.h>
#endif
// clang-format on


namespace AMP {


/*
 * Routine to convert an integer to a string.
 */
std::string Utilities::intToString( int num, int min_width )
{
    int tmp_width = ( min_width > 0 ? min_width : 1 );
    std::ostringstream os;
    if ( num < 0 ) {
        os << '-' << std::setw( tmp_width - 1 ) << std::setfill( '0' ) << -num;
    } else {
        os << std::setw( tmp_width ) << std::setfill( '0' ) << num;
    }
    os << std::flush;

    return ( os.str() ); // returns the string form of the stringstream object
}

std::string Utilities::nodeToString( int num ) { return intToString( num, 5 ); }

std::string Utilities::processorToString( int num ) { return intToString( num, 5 ); }


std::string Utilities::patchToString( int num ) { return intToString( num, 4 ); }


std::string Utilities::levelToString( int num ) { return intToString( num, 4 ); }


std::string Utilities::blockToString( int num ) { return intToString( num, 4 ); }


/*
 * Routine that calls abort and prints calling location to error stream.
 */
void Utilities::abort( const std::string &message, const std::string &filename, const int line )
{
    if ( AMP::AMPManager::use_MPI_Abort == true ) {
        // Get the call stack and memory usage
        long long unsigned int N_bytes = getMemoryUsage();
        auto stack = StackTrace::getCallStack();
        // Log the abort message
        if ( AMPManager::isInitialized() ) {
            AMP::pout << "Bytes used = " << N_bytes << std::endl;
            AMP::pout << "Stack Trace:\n";
            for ( auto &elem : stack )
                AMP::pout << "   " << elem.print() << std::endl;
            AMP::pout << std::endl;
            Logger::getInstance()->logAbort( message, filename, line );
        } else {
            std::cout << "Bytes used = " << N_bytes << std::endl;
            std::cout << "Stack Trace:\n";
            for ( auto &elem : stack )
                std::cout << "   " << elem.print() << std::endl;
            std::cout << std::endl;
            std::cerr << "Program abort called in file ``" << filename << "'' at line " << line << std::endl;
            std::cerr << "ERROR MESSAGE: " << std::endl << message.c_str() << std::endl;
            std::cerr << std::flush;
        }
        // Use MPI_abort (will terminate all processes)
        AMP_MPI comm = AMP_MPI( AMP_COMM_WORLD );
        comm.abort();
    } else {
        // Throw a standard exception (allows the use of try, catch)
        // std::stringstream  stream;
        // stream << message << std::endl << "  " << filename << ":  " << line;
        // std::cout << stream.str() << std::endl;
        throw std::logic_error( message );
    }
}


// Function to create a 32-bit hash key from a character array
unsigned int Utilities::hash_char( const char *name )
{
    AMP_INSIST( sizeof( unsigned int ) == 4, "Need unsigned 32-bit int" );
    unsigned int hash = 5381;
    unsigned char c;
    while ( ( c = *name++ ) ) {
        // hash = hash * 33 ^ c
        hash = ( ( hash << 5 ) + hash ) ^ c;
    }
    return hash;
}


/****************************************************************************
 *  Function to set an environemental variable                               *
 ****************************************************************************/
void Utilities::setenv( const char *name, const char *value )
{
#if defined( USE_LINUX ) || defined( USE_MAC )
    char env[100];
    sprintf( env, "%s=%s", name, value );
    bool pass = false;
    if ( value != nullptr )
        pass = ::setenv( name, value, 1 ) == 0;
    else
        pass = ::unsetenv( name ) == 0;
#elif defined( USE_WINDOWS )
    bool pass = SetEnvironmentVariable( name, value ) != 0;
#else
#error Unknown OS
#endif
    if ( !pass ) {
        char msg[100];
        if ( value != nullptr )
            sprintf( msg, "Error setting enviornmental variable: %s=%s\n", name, value );
        else
            sprintf( msg, "Error clearing enviornmental variable: %s\n", name );
        AMP_ERROR( msg );
    }
}


/****************************************************************************
 *  Function to get the memory usage                                         *
 *  Note: this function should be thread-safe                                *
 ****************************************************************************/
// clang-format off
#if defined( USE_MAC ) || defined( USE_LINUX )
    // Get the page size on mac or linux
    static size_t page_size = static_cast<size_t>( sysconf( _SC_PAGESIZE ) );
#endif
size_t Utilities::getSystemMemory()
{
    #if defined( USE_LINUX )
        static long pages = sysconf( _SC_PHYS_PAGES );
        size_t N_bytes    = pages * page_size;
    #elif defined( USE_MAC )
        int mib[2]    = { CTL_HW, HW_MEMSIZE };
        u_int namelen = sizeof( mib ) / sizeof( mib[0] );
        uint64_t size;
        size_t len = sizeof( size );
        size_t N_bytes = 0;
        if ( sysctl( mib, namelen, &size, &len, nullptr, 0 ) == 0 )
            N_bytes = size;
    #elif defined( USE_WINDOWS )
        MEMORYSTATUSEX status;
        status.dwLength = sizeof( status );
        GlobalMemoryStatusEx( &status );
        size_t N_bytes = status.ullTotalPhys;
    #else
        #error Unknown OS
    #endif
    return N_bytes;
}
size_t Utilities::getMemoryUsage()
{
    #ifdef USE_TIMER
        size_t N_bytes = MemoryApp::getTotalMemoryUsage();
    #else
        #if defined( USE_LINUX )
            struct mallinfo meminfo = mallinfo();
            size_t size_hblkhd      = static_cast<unsigned int>( meminfo.hblkhd );
            size_t size_uordblks    = static_cast<unsigned int>( meminfo.uordblks );
            size_t N_bytes          = size_hblkhd + size_uordblks;
        #elif defined( USE_MAC )
            struct task_basic_info t_info;
            mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
            if ( KERN_SUCCESS !=
                 task_info( mach_task_self(), TASK_BASIC_INFO, (task_info_t) &t_info, &t_info_count ) ) {
                return 0;
            }
            size_t N_bytes = t_info.virtual_size;
        #elif defined( USE_WINDOWS )
            PROCESS_MEMORY_COUNTERS memCounter;
            GetProcessMemoryInfo( GetCurrentProcess(), &memCounter, sizeof( memCounter ) );
            size_t N_bytes = memCounter.WorkingSetSize;
        #else
            #error Unknown OS
        #endif
    #endif
    return N_bytes;
}
// clang-format on


/****************************************************************************
 *  Functions to get the time and timer resolution                           *
 ****************************************************************************/
#if defined( USE_WINDOWS )
double Utilities::time()
{
    LARGE_INTEGER end, f;
    QueryPerformanceFrequency( &f );
    QueryPerformanceCounter( &end );
    double time = ( (double) end.QuadPart ) / ( (double) f.QuadPart );
    return time;
}
double Utilities::tick()
{
    LARGE_INTEGER f;
    QueryPerformanceFrequency( &f );
    double resolution = ( (double) 1.0 ) / ( (double) f.QuadPart );
    return resolution;
}
void Utilities::sleepMs( unsigned int N ) { Sleep( N ); }
#elif defined( USE_LINUX ) || defined( USE_MAC )
double Utilities::time()
{
    timeval current_time;
    gettimeofday( &current_time, nullptr );
    double time = ( (double) current_time.tv_sec ) + 1e-6 * ( (double) current_time.tv_usec );
    return time;
}
double Utilities::tick()
{
    timeval start, end;
    gettimeofday( &start, nullptr );
    gettimeofday( &end, nullptr );
    while ( end.tv_sec == start.tv_sec && end.tv_usec == start.tv_usec )
        gettimeofday( &end, nullptr );
    double resolution = ( (double) ( end.tv_sec - start.tv_sec ) ) +
                        1e-6 * ( (double) ( end.tv_usec - start.tv_usec ) );
    return resolution;
}
void Utilities::sleepMs( unsigned int N ) { usleep( N * 1000 ); }
#else
#error Unknown OS
#endif


/****************************************************************************
 *  Filesystem utilities                                                     *
 ****************************************************************************/
bool Utilities::fileExists( const std::string &filename )
{
    std::ifstream ifile( filename.c_str() );
    return ifile.good();
}

void Utilities::renameFile( const std::string &old_filename, const std::string &new_filename )
{
    AMP_ASSERT( !old_filename.empty() );
    AMP_ASSERT( !new_filename.empty() );
    rename( old_filename.c_str(), new_filename.c_str() );
}

void Utilities::deleteFile( const std::string &filename )
{
    AMP_ASSERT( !filename.empty() );
    if ( fileExists( filename ) ) {
        int error = remove( filename.c_str() );
        AMP_INSIST( error == 0, "Error deleting file" );
    }
}
void Utilities::recursiveMkdir( const std::string &path, mode_t mode, bool only_node_zero_creates )
{
    AMP_MPI comm = AMP_MPI( AMP_COMM_WORLD );
    if ( ( !only_node_zero_creates ) || ( comm.getRank() == 0 ) ) {
        auto length    = (int) path.length();
        auto *path_buf = new char[length + 1];
        sprintf( path_buf, "%s", path.c_str() );
        struct stat status;
        int pos = length - 1;
        // find part of path that has not yet been created
        while ( ( stat( path_buf, &status ) != 0 ) && ( pos >= 0 ) ) {
            // slide backwards in string until next slash found
            bool slash_found = false;
            while ( ( !slash_found ) && ( pos >= 0 ) ) {
                if ( path_buf[pos] == '/' || path_buf[pos] == 92 ) {
                    slash_found = true;
                    if ( pos >= 0 )
                        path_buf[pos] = '\0';
                } else
                    pos--;
            }
        }
        // if there is a part of the path that already exists make sure it is really a directory
        if ( pos >= 0 ) {
            if ( !S_ISDIR( status.st_mode ) ) {
                AMP_ERROR( "Error in Utilities::recursiveMkdir...\n"
                           << "    Cannot create directories in path = " << path
                           << "\n    because some intermediate item in path exists and"
                           << "is NOT a directory" << std::endl );
            }
        }
        // make all directories that do not already exist
        if ( pos < 0 ) {
            if ( mkdir( path_buf, mode ) != 0 ) {
                AMP_ERROR( "Error in Utilities::recursiveMkdir...\n"
                           << "    Cannot create directory  = " << path_buf << std::endl );
            }
            pos = 0;
        }
        // make rest of directories
        do {
            // slide forward in string until next '\0' found
            bool null_found = false;
            while ( ( !null_found ) && ( pos < length ) ) {
                if ( path_buf[pos] == '\0' ) {
                    null_found    = true;
                    path_buf[pos] = '/';
                }
                pos++;
            }
            // make directory if not at end of path
            if ( pos < length ) {
                if ( mkdir( path_buf, mode ) != 0 ) {
                    AMP_ERROR( "Error in Utilities::recursiveMkdir...\n"
                               << "    Cannot create directory  = " << path_buf << std::endl );
                }
            }
        } while ( pos < length );
        delete[] path_buf;
    }
    // Make sure all processors wait until node zero creates the directory structure.
    if ( only_node_zero_creates )
        comm.barrier();
}
 
 
/****************************************************************************
 *  Print AMP Banner                                                         *
 ****************************************************************************/
void Utilities::printBanner()
{
    std::ostringstream banner;
    banner << std::endl;
    banner << "            _____                    _____                    _____" << std::endl;
    banner << R"(           /\    \                  /\    \                  /\    \ )"
           << std::endl;
    banner << R"(          /::\    \                /::\____\                /::\    \)"
           << std::endl;
    banner << R"(         /::::\    \              /::::|   |               /::::\    \)"
           << std::endl;
    banner << R"(        /::::::\    \            /:::::|   |              /::::::\    \)"
           << std::endl;
    banner << R"(       /:::/\:::\    \          /::::::|   |             /:::/\:::\    \)"
           << std::endl;
    banner << R"(      /:::/__\:::\    \        /:::/|::|   |            /:::/__\:::\    \)"
           << std::endl;
    banner << R"(     /::::\   \:::\    \      /:::/ |::|   |           /::::\   \:::\    \)"
           << std::endl;
    banner << R"(    /::::::\   \:::\    \    /:::/  |::|___|______    /::::::\   \:::\    \)"
           << std::endl;
    banner << "   /:::/\\:::\\   \\:::\\    \\  /:::/   |::::::::\\    \\  /:::/\\:::\\   "
              "\\:::\\____\\"
           << std::endl;
    banner << R"(  /:::/  \:::\   \:::\____\/:::/    |:::::::::\____\/:::/  \:::\   \:::|    |)"
           << std::endl;
    banner << R"(  \::/    \:::\  /:::/    /\::/    / ~~~~~/:::/    /\::/    \:::\  /:::|____|)"
           << std::endl;
    banner << R"(   \/____/ \:::\/:::/    /  \/____/      /:::/    /  \/_____/\:::\/:::/    /)"
           << std::endl;
    banner << "            \\::::::/    /               /:::/    /            \\::::::/    /"
           << std::endl;
    banner << "             \\::::/    /               /:::/    /              \\::::/    /"
           << std::endl;
    banner << "             /:::/    /               /:::/    /                \\::/____/"
           << std::endl;
    banner << "            /:::/    /               /:::/    /" << std::endl;
    banner << "           /:::/    /               /:::/    /" << std::endl;
    banner << "          /:::/    /               /:::/    /" << std::endl;
    banner << "          \\::/    /                \\::/    /" << std::endl;
    banner << "           \\/____/                  \\/____/" << std::endl;
    banner << std::endl << std::endl;

    AMP::pout << banner.str();
}

// Factor a number into it's prime factors
std::vector<int> Utilities::factor( size_t number )
{
    if ( number <= 3 )
        return std::vector<int>( 1, (int) number );
    size_t i, n, n_max;
    bool factor_found;
    // Compute the maximum number of factors
    int N_primes_max = 1;
    n                = number;
    while ( n >>= 1 )
        ++N_primes_max;
    // Initialize n, factors
    n = number;
    std::vector<int> factors;
    factors.reserve( N_primes_max );
    while ( true ) {
        // Check if n is a trivial prime number
        if ( n == 2 || n == 3 || n == 5 ) {
            factors.push_back( (int) n );
            break;
        }
        // Check if n is divisible by 2
        if ( n % 2 == 0 ) {
            factors.push_back( 2 );
            n /= 2;
            continue;
        }
        // Check each odd number until a factor is reached
        n_max        = (size_t) floor( sqrt( (double) n ) );
        factor_found = false;
        for ( i = 3; i <= n_max; i += 2 ) {
            if ( n % i == 0 ) {
                factors.push_back( i );
                n /= i;
                factor_found = true;
                break;
            }
        }
        if ( factor_found )
            continue;
        // No factors were found, the number must be prime
        factors.push_back( (int) n );
        break;
    }
    // Sort the factors
    AMP::Utilities::quicksort( factors );
    return factors;
}


// Function to perform linear interpolation
double Utilities::linear( const std::vector<double> &x, const std::vector<double> &f, double xi )
{
    size_t Nx = x.size();
    AMP_ASSERT( Nx > 1 );
    AMP_ASSERT( f.size() == Nx );
    size_t i = AMP::Utilities::findfirst( x, xi );
    if ( i == 0 ) {
        i = 1;
    }
    if ( i == x.size() ) {
        i = x.size() - 1;
    }
    double dx = ( xi - x[i - 1] ) / ( x[i] - x[i - 1] );
    return dx * f[i] + ( 1.0 - dx ) * f[i - 1];
}


// Function to perform bi-linear interpolation
double Utilities::bilinear( const std::vector<double> &x,
                            const std::vector<double> &y,
                            const std::vector<double> &f,
                            double xi,
                            double yi )
{
    size_t Nx = x.size();
    size_t Ny = y.size();
    AMP_ASSERT( Nx > 1 && Ny > 1 );
    AMP_ASSERT( f.size() == Nx * Ny );
    size_t i = AMP::Utilities::findfirst( x, xi );
    size_t j = AMP::Utilities::findfirst( y, yi );
    if ( i == 0 ) {
        i = 1;
    }
    if ( j == 0 ) {
        j = 1;
    }
    if ( i == x.size() ) {
        i = x.size() - 1;
    }
    if ( j == y.size() ) {
        j = y.size() - 1;
    }
    double dx  = ( xi - x[i - 1] ) / ( x[i] - x[i - 1] );
    double dy  = ( yi - y[j - 1] ) / ( y[j] - y[j - 1] );
    double f1  = f[i - 1 + ( j - 1 ) * Nx];
    double f2  = f[i + ( j - 1 ) * Nx];
    double f3  = f[i - 1 + j * Nx];
    double f4  = f[i + j * Nx];
    double dx2 = 1.0 - dx;
    double dy2 = 1.0 - dy;
    return ( dx * f2 + dx2 * f1 ) * dy2 + ( dx * f4 + dx2 * f3 ) * dy;
}


// Function to perform tri-linear interpolation
double Utilities::trilinear( const std::vector<double> &x,
                             const std::vector<double> &y,
                             const std::vector<double> &z,
                             const std::vector<double> &f,
                             double xi,
                             double yi,
                             double zi )
{
    size_t Nx = x.size();
    size_t Ny = y.size();
    size_t Nz = z.size();
    AMP_ASSERT( Nx > 1 && Ny > 1 && Nz > 1 );
    AMP_ASSERT( f.size() == Nx * Ny * Nz );
    size_t i = AMP::Utilities::findfirst( x, xi );
    size_t j = AMP::Utilities::findfirst( y, yi );
    size_t k = AMP::Utilities::findfirst( z, zi );
    if ( i == 0 ) {
        i = 1;
    }
    if ( j == 0 ) {
        j = 1;
    }
    if ( k == 0 ) {
        k = 1;
    }
    if ( i == x.size() ) {
        i = x.size() - 1;
    }
    if ( j == y.size() ) {
        j = y.size() - 1;
    }
    if ( k == z.size() ) {
        k = z.size() - 1;
    }
    double dx  = ( xi - x[i - 1] ) / ( x[i] - x[i - 1] );
    double dy  = ( yi - y[j - 1] ) / ( y[j] - y[j - 1] );
    double dz  = ( zi - z[k - 1] ) / ( z[k] - z[k - 1] );
    double f1  = f[i - 1 + ( j - 1 ) * Nx + ( k - 1 ) * Nx * Ny];
    double f2  = f[i + ( j - 1 ) * Nx + ( k - 1 ) * Nx * Ny];
    double f3  = f[i - 1 + j * Nx + ( k - 1 ) * Nx * Ny];
    double f4  = f[i + j * Nx + ( k - 1 ) * Nx * Ny];
    double f5  = f[i - 1 + ( j - 1 ) * Nx + k * Nx * Ny];
    double f6  = f[i + ( j - 1 ) * Nx + k * Nx * Ny];
    double f7  = f[i - 1 + j * Nx + k * Nx * Ny];
    double f8  = f[i + j * Nx + k * Nx * Ny];
    double dx2 = 1.0 - dx;
    double dy2 = 1.0 - dy;
    double dz2 = 1.0 - dz;
    double h0  = ( dx * f2 + dx2 * f1 ) * dy2 + ( dx * f4 + dx2 * f3 ) * dy;
    double h1  = ( dx * f6 + dx2 * f5 ) * dy2 + ( dx * f8 + dx2 * f7 ) * dy;
    return h0 * dz2 + h1 * dz;
}


// Dummy function to prevent compiler from optimizing away variable
void Utilities::nullUse( void *data ) { NULL_USE( data ); }


} // namespace AMP
