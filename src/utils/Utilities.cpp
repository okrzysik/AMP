#define NOMINMAX
#include "AMP/utils/Utilities.hpp"
#include "AMP/AMP_TPLs.h"
#include "AMP/IO/PIO.h"
#include "AMP/mesh/MeshID.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"

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

#ifdef AMP_USE_TIMER
    #include "MemoryApp.h"
#endif

#if defined( AMP_USE_CUDA ) || defined( USE_CUDA )
    #include "AMP/utils/cuda/helper_cuda.h"
#endif
#if defined( AMP_USE_HIP ) || defined( USE_HIP )
    #include "AMP/utils/hip/helper_hip.h"
#endif

// Include system dependent headers
// clang-format off
#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 ) || defined( _MSC_VER )
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
#if defined( __APPLE__ )
    #include <mach/mach.h>
    #include <sys/sysctl.h>
    #include <sys/types.h>
#endif
// clang-format on


namespace AMP::Utilities {


// Mutex for Utility functions
static std::mutex Utilities_mutex;


/*
 * Routine to convert an integer to a string.
 */
std::string intToString( int num, int min_width )
{
    int tmp_width = ( min_width > 0 ? min_width : 1 );
    std::ostringstream os;
    if ( num < 0 ) {
        os << '-' << std::setw( tmp_width - 1 ) << std::setfill( '0' ) << -num;
    } else {
        os << std::setw( tmp_width ) << std::setfill( '0' ) << num;
    }
    os << std::flush;
    return os.str();
}
std::string nodeToString( int num ) { return intToString( num, 5 ); }
std::string processorToString( int num ) { return intToString( num, 5 ); }
std::string patchToString( int num ) { return intToString( num, 4 ); }
std::string levelToString( int num ) { return intToString( num, 4 ); }
std::string blockToString( int num ) { return intToString( num, 4 ); }


/****************************************************************************
 *  Basic checks                                                             *
 ****************************************************************************/
static_assert( getOS() != OS::Unknown );


/****************************************************************************
 *  Get pointer location                                                     *
 ****************************************************************************/
MemoryType getMemoryType( const void *ptr )
{
    NULL_USE( ptr );
    auto type = MemoryType::unregistered;
#if defined( AMP_USE_CUDA ) || defined( USE_CUDA )
    type = getCudaMemoryType( ptr );
    if ( type != MemoryType::unregistered )
        return type;
#endif
#if defined( AMP_USE_HIP ) || defined( USE_HIP )
    type = getHipMemoryType( ptr );
    if ( type != MemoryType::unregistered )
        return type;
#endif
    return type;
}
std::string getString( MemoryType type )
{
    if ( type == MemoryType::unregistered )
        return "unregistered";
    else if ( type == MemoryType::host )
        return "host";
    else if ( type == MemoryType::device )
        return "device";
    else if ( type == MemoryType::managed )
        return "managed";
    else
        AMP_ERROR( "Unknown pointer type" );
}


/****************************************************************************
 *  Function to set an environemental variable                               *
 ****************************************************************************/
void setenv( const char *name, const char *value )
{
    Utilities_mutex.lock();
#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 ) || \
    defined( _MSC_VER )
    bool pass = SetEnvironmentVariable( name, value ) != 0;
#else
    bool pass = false;
    if ( value == nullptr )
        pass = ::unsetenv( name ) == 0;
    else
        pass = ::setenv( name, value, 1 ) == 0;
#endif
    Utilities_mutex.unlock();
    if ( !pass ) {
        std::string msg;
        if ( value != nullptr )
            msg = stringf( "Error setting enviornmental variable: %s=%s\n", name, value );
        else
            msg = stringf( "Error clearing enviornmental variable: %s\n", name );
        AMP_ERROR( msg );
    }
}
std::string getenv( const char *name )
{
    std::string var;
    Utilities_mutex.lock();
    auto tmp = std::getenv( name );
    if ( tmp )
        var = std::string( tmp );
    Utilities_mutex.unlock();
    return var;
}


/****************************************************************************
 *  Print AMP Banner                                                         *
 ****************************************************************************/
// clang-format off
void printBanner()
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


/****************************************************************************
 *  Prime number functions                                                   *
 ****************************************************************************/
std::vector<int> factor( uint64_t n )
{
    // Handle trival case
    if ( n <= 3 )
        return { static_cast<int>( n ) };
    // Initialize factors
    size_t N = 0;
    int factors[64];
    // Remove all factors of 2
    while ( ( n & 0x01 ) == 0 ) {
        factors[N++] = 2;
        n >>= 1;
    }
    // Remove all factors of 3
    while ( n % 3 == 0 ) {
        factors[N++] = 3;
        n /= 3;
    }
    // Use brute force to find remaining factors
    uint64_t f = 5;
    while ( true ) {
        // Determine the largest number we need to check
        auto f_max = static_cast<uint64_t>( floor( 1.000000000000001 * sqrt( n ) ) );
        // Search all remaining numbers (note  we skip every 3rd odd number)
        bool found = false;
        for ( ; f <= f_max && !found; f += 6 ) {
            while ( n % f == 0 ) {
                factors[N++] = f;
                n /= f;
                found = true;
            }
            while ( n % ( f + 2 ) == 0 ) {
                factors[N++] = f + 2;
                n /= f + 2;
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
bool isPrime( uint64_t n )
{
    if ( n <= 3 )
        return true;
    if ( ( n & 0x01 ) == 0 || n % 3 == 0 )
        return false;
    // Determine the largest number we need to check
    auto f_max = static_cast<uint64_t>( floor( 1.000000000000001 * sqrt( n ) ) );
    // Check if the number is prime
    for ( uint64_t f = 5; f <= f_max; f += 6 ) {
        if ( ( n % f == 0 ) || ( n % ( f + 2 ) == 0 ) )
            return false;
    }
    return true;
}
std::vector<uint64_t> primes( uint64_t n )
{
    // Handle special cases
    if ( n < 2 )
        return { 1u };
    if ( n == 2 )
        return { 2u };
    // Create our bit array
    uint64_t n2 = ( n + 1 ) / 2;
    double tmp  = 1.000000000000001 * sqrt( static_cast<double>( n ) );
    uint64_t ub = static_cast<uint64_t>( tmp ) >> 1;
    auto N      = ( n2 + 63 ) / 64;
    auto p      = new uint64_t[N];
    memset( p, 0xFF, sizeof( uint64_t ) * N );
    // Helper functions to get/set the bits
    auto get = [p]( size_t i ) {
        size_t i1 = i >> 6;
        size_t i2 = i & 0x3F;
        return ( p[i1] & ( 1UL << i2 ) ) != 0;
    };
    auto unset = [p]( size_t i ) {
        size_t i1 = i >> 6;
        size_t i2 = i & 0x3F;
        p[i1] &= ~( 1UL << i2 );
    };
    // Set all non-prime values to false
    static_assert( sizeof( unsigned long ) == sizeof( uint64_t ) );
    for ( uint64_t k = 1; k <= ub; k++ ) {
        if ( get( k ) ) {
            uint64_t k2 = 2 * k + 1;
            for ( size_t j = 2 * k * ( k + 1 ); j < n2; j += k2 )
                unset( j );
        }
    }
    // Store the prime numbers (note: this takes longer than computing them)
    auto M = static_cast<size_t>( n / log2( n ) );
    M      = 1UL << static_cast<int>( round( log2( M ) ) );
    std::vector<uint64_t> p2;
    p2.reserve( M );
    p2.push_back( 2 );
    for ( size_t i = 1; i < n2; i++ ) {
        if ( get( i ) )
            p2.push_back( 2 * i + 1 );
    }
    return p2;
}


// Function to perform linear interpolation
double linear( const std::vector<double> &x, const std::vector<double> &f, double xi )
{
    size_t Nx = x.size();
    AMP_ASSERT( Nx > 1 );
    AMP_ASSERT( f.size() == Nx );
    size_t i = findfirst( x, xi );
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
double bilinear( const std::vector<double> &x,
                 const std::vector<double> &y,
                 const std::vector<double> &f,
                 double xi,
                 double yi )
{
    size_t Nx = x.size();
    size_t Ny = y.size();
    AMP_ASSERT( Nx > 1 && Ny > 1 );
    AMP_ASSERT( f.size() == Nx * Ny );
    size_t i = findfirst( x, xi );
    size_t j = findfirst( y, yi );
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
double trilinear( const std::vector<double> &x,
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
    size_t i = findfirst( x, xi );
    size_t j = findfirst( y, yi );
    size_t k = findfirst( z, zi );
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
void nullUse( void *data ) { NULL_USE( data ); }


// Function to demangle a string (e.g. from typeid)
#ifdef __GNUC__
    #define USE_ABI
    #include <cxxabi.h>
#endif
std::string demangle( const std::string &name )
{
    std::string out;
#ifdef __GNUC__
    int status;
    char *demangled = abi::__cxa_demangle( name.data(), nullptr, nullptr, &status );
    if ( status == 0 && demangled != nullptr )
        out = demangled;
    free( demangled );
#endif
    if ( out.empty() )
        out = name;
    return out;
}


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
void printDatabase( const Database &db, std::ostream &os, const std::string &indent )
{
    db.print( os, indent );
}


} // namespace AMP::Utilities


/************************************************************************
 * Explicit instantiations                                               *
 ************************************************************************/
#define INSTANTIATE( TYPE )                                                            \
    template std::string AMP::Utilities::to_string<TYPE>( std::vector<TYPE> const & ); \
    template void AMP::Utilities::quicksort<TYPE>( size_t, TYPE * );                   \
    template void AMP::Utilities::quicksort<TYPE, TYPE>( size_t, TYPE *, TYPE * );     \
    template void AMP::Utilities::unique<TYPE>( std::vector<TYPE> & );                 \
    template void AMP::Utilities::unique<TYPE>(                                        \
        std::vector<TYPE> &, std::vector<size_t> &, std::vector<size_t> & );           \
    template size_t AMP::Utilities::findfirst<TYPE>( size_t, const TYPE *, const TYPE & )
INSTANTIATE( int8_t );
INSTANTIATE( int16_t );
INSTANTIATE( int32_t );
INSTANTIATE( int64_t );
INSTANTIATE( uint8_t );
INSTANTIATE( uint16_t );
INSTANTIATE( uint32_t );
INSTANTIATE( uint64_t );
INSTANTIATE( float );
INSTANTIATE( double );
INSTANTIATE( long double );
