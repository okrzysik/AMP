#define NOMINMAX
#include "AMP/utils/Utilities.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
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
#include <mutex>
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
#elif defined( __linux ) || defined( __linux__ ) || defined( __unix ) || defined( __posix )
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


// Mutex for Utility functions
static std::mutex Utilities_mutex;


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
    return ( os.str() );
}
std::string Utilities::nodeToString( int num ) { return intToString( num, 5 ); }
std::string Utilities::processorToString( int num ) { return intToString( num, 5 ); }
std::string Utilities::patchToString( int num ) { return intToString( num, 4 ); }
std::string Utilities::levelToString( int num ) { return intToString( num, 4 ); }
std::string Utilities::blockToString( int num ) { return intToString( num, 4 ); }


/****************************************************************************
 *  Function to set an environemental variable                               *
 ****************************************************************************/
void Utilities::setenv( const char *name, const char *value )
{
    Utilities_mutex.lock();
#if defined( USE_LINUX ) || defined( USE_MAC )
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
    Utilities_mutex.unlock();
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
 *  Filesystem utilities                                                     *
 ****************************************************************************/
bool Utilities::fileExists( const std::string &filename )
{
    std::ifstream ifile( filename.c_str() );
    return ifile.good();
}

std::string Utilities::path( const std::string &filename )
{
    size_t pos = 0;
    if ( filename.find_last_of( 47 ) != std::string::npos )
        pos = filename.find_last_of( 47 );
    if ( filename.find_last_of( 92 ) != std::string::npos )
        pos = std::max( pos, filename.find_last_of( 92 ) );
    return filename.substr( 0, pos );
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

std::string Utilities::getSuffix( const std::string &filename )
{
    size_t pos = filename.rfind( '.' );
    if ( pos == std::string::npos )
        return std::string();
    std::string suffix = filename.substr( pos + 1 );
    for ( auto &c : suffix )
        c = std::tolower( c );
    return suffix;
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
// clang-format off
void Utilities::printBanner()
{
    constexpr char banner[] =
        R"(            _____                    _____                    _____)" "\n"
        R"(           /\    \                  /\    \                  /\    \)" "\n"
        R"(          /::\    \                /::\____\                /::\    \)" "\n"
        R"(         /::::\    \              /::::|   |               /::::\    \)" "\n"
        R"(        /::::::\    \            /:::::|   |              /::::::\    \)" "\n"
        R"(       /:::/\:::\    \          /::::::|   |             /:::/\:::\    \)" "\n"
        R"(      /:::/__\:::\    \        /:::/|::|   |            /:::/__\:::\    \)" "\n"
        R"(     /::::\   \:::\    \      /:::/ |::|   |           /::::\   \:::\    \)" "\n"
        R"(    /::::::\   \:::\    \    /:::/  |::|___|______    /::::::\   \:::\    \)" "\n"
        R"(   /:::/\:::\   \:::\    \  /:::/   |::::::::\    \  /:::/\:::\   \:::\____\)" "\n"
        R"(  /:::/  \:::\   \:::\____\/:::/    |:::::::::\____\/:::/  \:::\   \:::|    |)" "\n"
        R"(  \::/    \:::\  /:::/    /\::/    / ~~~~~/:::/    /\::/    \:::\  /:::|____|)" "\n"
        R"(   \/____/ \:::\/:::/    /  \/____/      /:::/    /  \/_____/\:::\/:::/    /)" "\n"
        R"(            \::::::/    /               /:::/    /            \::::::/    /)" "\n"
        R"(             \::::/    /               /:::/    /              \::::/    /)" "\n"
        R"(             /:::/    /               /:::/    /                \::/____/)" "\n"
        R"(            /:::/    /               /:::/    /)" "\n"
        R"(           /:::/    /               /:::/    /)" "\n"
        R"(          /:::/    /               /:::/    /)" "\n"
        R"(          \::/    /                \::/    /)" "\n"
        R"(           \/____/                  \/____/)" "\n";
    AMP::pout << "\n" << banner << "\n" << std::endl;
}
// clang-format on

// Factor a number into it's prime factors
std::vector<int> Utilities::factor( uint64_t number )
{
    uint64_t n = number;
    // Handle trival case
    if ( n <= 3 )
        return std::vector<int>( 1, (int) number );
    // Initialize factors
    size_t N = 0;
    int factors[64];
    // Remove all factors of 2
    while ( ( n & 0x01 ) == 0 ) {
        factors[N++] = 2;
        n >>= 1;
        continue;
    }
    // Use brute force to find remaining factors
    uint64_t f = 3;
    while ( true ) {
        auto f_max = static_cast<uint64_t>( floor( sqrt( n ) ) );
        bool found = false;
        for ( ; f <= f_max && !found; f += 2 ) {
            while ( n % f == 0 ) {
                factors[N++] = f;
                n /= f;
                found = true;
            }
        }
        if ( !found ) {
            factors[N++] = n;
            break;
        }
    }
    return std::vector<int>( factors, factors + N );
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


// Print a database to an output stream
template<class TYPE>
static void printVar( const std::string &name,
                      const std::vector<TYPE> &data,
                      std::ostream &os,
                      const std::string &indent )
{
    os << indent << name << " = ";
    if ( !data.empty() ) {
        os << data[0];
        for ( size_t i = 1; i < data.size(); i++ )
            os << ", " << data[i];
    }
    os << std::endl;
}
void Utilities::printDatabase( Database &db, std::ostream &os, const std::string &indent )
{
    for ( const auto &name : db.getAllKeys() ) {
        if ( db.isDatabase( name ) ) {
            os << indent << name << "{\n";
            printDatabase( *db.getDatabase( name ), os, indent + "   " );
            os << indent << "}\n";
        } else if ( db.isType<bool>( name ) ) {
            auto data = db.getVector<bool>( name );
            os << indent << name << " = ";
            if ( !data.empty() ) {
                os << ( data[0] ? "TRUE" : "FALSE" );
                for ( size_t i = 1; i < data.size(); i++ )
                    os << ", " << ( data[i] ? "TRUE" : "FALSE" );
            }
            os << std::endl;
        } else if ( db.isType<double>( name ) ) {
            auto data = db.getVector<double>( name );
            printVar( name, data, os, indent );
        } else if ( db.isType<int>( name ) ) {
            auto data = db.getVector<int>( name );
            printVar( name, data, os, indent );
        } else if ( db.isString( name ) ) {
            auto data = db.getVector<std::string>( name );
            printVar( name, data, os, indent );
        } else {
            AMP_ERROR( "Unknown type for field: " + name );
        }
    }
}

} // namespace AMP
