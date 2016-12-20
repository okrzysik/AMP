#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>
#include <memory>

#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/StackTrace.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/enable_shared_from_this.h"
#include "utils/shared_ptr.h"


// Detect the OS (defines which tests we allow to fail)
#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 ) || \
    defined( _MSC_VER )
#define USE_WINDOWS
#elif defined( __APPLE__ )
#define USE_MAC
#elif defined( __linux ) || defined( __unix ) || defined( __posix )
#define USE_LINUX
#else
#error Unknown OS
#endif


using namespace AMP;


// Subract two size_t numbers, returning the absolute value
size_t abs_diff( size_t a, size_t b ) { return ( a >= b ) ? a - b : b - a; }


// This checks approx_equal
template <class T>
void testApproxEqualInt( UnitTest *ut )
{
    std::string type_name( typeid( T ).name() );
    if ( Utilities::approx_equal<T>( 100000, 100000 ) &&
         Utilities::approx_equal_abs<T>( 100000, 100000 ) &&
         !Utilities::approx_equal<T>( 100000, 100001 ) &&
         !Utilities::approx_equal_abs<T>( 100000, 100001 ) )
        ut->passes( "Integer (" + type_name + ") passes simple check." );
    else
        ut->failure( "Integer (" + type_name + ") passes simple check." );

    if ( Utilities::approx_equal_abs<T>( 100001, 100000, 1 ) &&
         !Utilities::approx_equal_abs<T>( 100002, 100000, 1 ) )
        ut->passes( "Integer (" + type_name + ") passes close simple check." );
    else
        ut->failure( "Integer (" + type_name + ") passes close simple check." );
}
template <class T>
void testApproxEqual( UnitTest *ut )
{
    std::string type_name( typeid( T ).name() );

    T mine      = 1.0;
    T close_rel = mine * static_cast<T>( 1.0 + pow( std::numeric_limits<T>::epsilon(), (T) 0.8 ) );
    T wrong_rel = mine * static_cast<T>( 1.0 + pow( std::numeric_limits<T>::epsilon(), (T) 0.7 ) );
    T close_abs = mine + pow( std::numeric_limits<T>::epsilon(), (T) 0.8 );
    T wrong_abs = mine + pow( std::numeric_limits<T>::epsilon(), (T) 0.7 );
    if ( Utilities::approx_equal( mine, close_rel ) &&
         Utilities::approx_equal_abs( mine, close_abs ) &&
         !Utilities::approx_equal( mine, wrong_rel ) &&
         !Utilities::approx_equal_abs( mine, wrong_abs ) )
        ut->passes( type_name + " passes simple check near 1" );
    else
        ut->failure( type_name + " passes simple check near 1" );

    mine      = static_cast<T>( 1e-6 );
    close_rel = mine * static_cast<T>( 1.0 + pow( std::numeric_limits<T>::epsilon(), (T) 0.8 ) );
    wrong_rel = mine * static_cast<T>( 1.0 + pow( std::numeric_limits<T>::epsilon(), (T) 0.7 ) );
    close_abs = mine + pow( std::numeric_limits<T>::epsilon(), (T) 0.8 );
    wrong_abs = mine + pow( std::numeric_limits<T>::epsilon(), (T) 0.7 );
    if ( Utilities::approx_equal( mine, close_rel ) &&
         Utilities::approx_equal_abs( mine, close_abs ) &&
         !Utilities::approx_equal( mine, wrong_rel ) &&
         !Utilities::approx_equal_abs( mine, wrong_abs ) )
        ut->passes( type_name + " passes simple check near 1e-6" );
    else
        ut->failure( type_name + " passes simple check near 1e-6" );

    mine      = static_cast<T>( -1e-32 );
    close_rel = mine * static_cast<T>( 1.0 + pow( std::numeric_limits<T>::epsilon(), (T) 0.8 ) );
    wrong_rel = mine * static_cast<T>( 1.0 + pow( std::numeric_limits<T>::epsilon(), (T) 0.7 ) );
    close_abs = mine + pow( std::numeric_limits<T>::epsilon(), (T) 0.8 );
    wrong_abs = mine + pow( std::numeric_limits<T>::epsilon(), (T) 0.7 );
    if ( Utilities::approx_equal( mine, close_rel ) &&
         Utilities::approx_equal_abs( mine, close_abs ) &&
         !Utilities::approx_equal( mine, wrong_rel ) &&
         !Utilities::approx_equal_abs( mine, wrong_abs ) )
        ut->passes( type_name + " passes simple check near -1e-32" );
    else
        ut->failure( type_name + " passes simple check near -1e-32" );
}


// Function to return the call stack
std::vector<StackTrace::stack_info> get_call_stack()
{
    return StackTrace::getCallStack();
}


// Function to test the interpolants
void test_interp( UnitTest *ut )
{
    const double a  = 1.0;
    const double bx = 1.0;
    const double by = -1.0;
    const double bz = 0.5;
    int Nx          = 20;
    int Ny          = 10;
    int Nz          = 5;
    std::vector<double> x( Nx, 0.0 );
    std::vector<double> y( Ny, 0.0 );
    std::vector<double> z( Nz, 0.0 );
    std::vector<double> f1( Nx, 0.0 );
    std::vector<double> f2( Nx * Ny, 0.0 );
    std::vector<double> f3( Nx * Ny * Nz, 0.0 );
    for ( int i = 0; i < Nx; i++ ) {
        x[i]  = ( (double) i ) / ( (double) ( Nx - 1 ) );
        f1[i] = a + bx * x[i];
        for ( int j = 0; j < Ny; j++ ) {
            y[j]           = ( (double) j ) / ( (double) ( Ny - 1 ) );
            f2[i + j * Nx] = a + bx * x[i] + by * y[j];
            for ( int k = 0; k < Nz; k++ ) {
                z[k]                         = ( (double) k ) / ( (double) ( Nz - 1 ) );
                f3[i + j * Nx + k * Nx * Ny] = a + bx * x[i] + by * y[j] + bz * z[k];
            }
        }
    }
    bool pass_linear    = true;
    bool pass_bilinear  = true;
    bool pass_trilinear = true;
    int Nix             = 100;
    int Niy             = 200;
    int Niz             = 50;
    for ( int i = 0; i < Nix; i++ ) {
        double xi = ( (double) i - 2 ) / ( (double) ( Nix - 5 ) );
        double fi = Utilities::linear( x, f1, xi );
        if ( !Utilities::approx_equal( fi, a + bx * xi, 1e-12 ) )
            pass_linear = false;
        for ( int j = 0; j < Niy; j++ ) {
            double yi = ( (double) j - 2 ) / ( (double) ( Niy - 5 ) );
            fi        = Utilities::bilinear( x, y, f2, xi, yi );
            if ( !Utilities::approx_equal( fi, a + bx * xi + by * yi, 1e-12 ) )
                pass_bilinear = false;
            for ( int k = 0; k < Niz; k++ ) {
                double zi = ( (double) k - 2 ) / ( (double) ( Niz - 5 ) );
                fi        = Utilities::trilinear( x, y, z, f3, xi, yi, zi );
                if ( !Utilities::approx_equal( fi, a + bx * xi + by * yi + bz * zi, 1e-12 ) )
                    pass_trilinear = false;
            }
        }
    }
    if ( pass_linear )
        ut->passes( "Linear interpolation" );
    else
        ut->failure( "Linear interpolation" );
    if ( pass_bilinear )
        ut->passes( "Bi-linear interpolation" );
    else
        ut->failure( "Bi-linear interpolation" );
    if ( pass_trilinear )
        ut->passes( "Tri-linear interpolation" );
    else
        ut->failure( "Tri-linear interpolation" );
}


// Test enable_shared_from_this
class dummy : public AMP::enable_shared_from_this<dummy>
{
public:
    shared_ptr<dummy> getPtr() { return shared_from_this(); }
};
static inline bool test_shared_from_this_pointer( const shared_ptr<dummy>& p1 )
{
    bool pass = p1.use_count() == 1;
    int count1, count2;
    auto p2   = p1->getPtr();
    count1    = p1.use_count();
    count2    = p2.use_count();
    pass      = pass && count1 == 2 && count2 == 2;
    shared_ptr<dummy> p3( p1 );
    count1  = p2.use_count();
    count2  = p3.use_count();
    pass    = pass && count1 == 3 && count2 == 3;
    p2.reset();
    count1  = p2.use_count();
    count2  = p3.use_count();
    pass    = pass && count1 == 0 && count2 == 2;
    auto p4 = p3->getPtr();
    count1  = p3.use_count();
    count2  = p4.use_count();
    pass    = pass && count1 == 3 && count2 == 3;
    shared_ptr<dummy> p5( p3.get(), []( void * ) {} );
    count1  = p3.use_count();
    count2  = p5.use_count();
    pass    = pass && count1 == 3 && count2 == 1;
    p5.reset();
    count1  = p3.use_count();
    count2  = p5.use_count();
    pass    = pass && p3.use_count() == 3 && p5.use_count() == 0;
    return pass;
}
void test_shared_from_this( UnitTest *ut )
{
    bool pass = true;
    try {
        shared_ptr<dummy> ptr( new dummy );
        pass = test_shared_from_this_pointer( ptr );
    } catch ( ... ) {
        pass = false;
    }
    if ( pass )
        ut->passes( "shared_from_this 1" );
    else
        ut->failure( "shared_from_this 1" );
    try {
        dummy *p1 = new dummy;
        auto p2   = p1->getPtr();
        pass = test_shared_from_this_pointer( p2 );
    } catch ( ... ) {
        pass = false;
    }
    if ( pass )
        ut->passes( "shared_from_this 2" );
    else
        ut->failure( "shared_from_this 2" );
}


/****************************************************************
* Run some basic utility tests                                  *
****************************************************************/
int main( int argc, char *argv[] )
{

    // Control the behavior of the startup
    AMPManagerProperties startup_properties;

    // Start AMP
    AMPManager::startup( argc, argv, startup_properties );
    int num_failed = 0;

    // Limit the scope of variables
    {
        AMP_MPI globalComm( AMP_COMM_WORLD );

        // Create the unit test
        UnitTest ut;

        // Print the banner
        Utilities::printBanner();

        // Test enable_shared_from_this
        test_shared_from_this( &ut );

        // Try converting an int to a string
        if ( Utilities::intToString( 37, 0 ) == "37" &&
             Utilities::intToString( 37, 3 ) == "037" )
            ut.passes( "Convert int to string" );
        else
            ut.failure( "Convert int to string" );

        // Test approx_equal
        testApproxEqualInt<int>( &ut );
        testApproxEqualInt<unsigned int>( &ut );
        testApproxEqualInt<size_t>( &ut );
        testApproxEqual<float>( &ut );
        testApproxEqual<double>( &ut );

        // Test interpolations
        test_interp( &ut );

        // Test quicksort performance
        size_t N = 10000;
        std::vector<int> data1( N );
        srand( static_cast<unsigned int>( time( nullptr ) ) );
        for ( size_t i         = 0; i < N; i++ )
            data1[i]           = rand();
        std::vector<int> data2 = data1;
        std::vector<int> data3 = data1;
        double t1              = Utilities::time();
        Utilities::quicksort( data1 );
        double t2 = Utilities::time();
        std::sort( data2.begin(), data2.end() );
        double t3 = Utilities::time();
        std::sort( &data3[0], &data3[0] + data3.size() );
        double t4 = Utilities::time();
        bool pass = true;
        for ( size_t i = 0; i < N; i++ ) {
            if ( data1[i] != data2[i] )
                pass = false;
        }
        if ( pass )
            ut.passes( "quicksort sorts correctly" );
        else
            ut.failure( "quicksort sorts correctly" );
        std::cout << "quicksort = " << t2 - t1 << ", std::sort = " << t3 - t2
                  << ", std::sort(2) = " << t4 - t3 << std::endl;

        // Test the hash key
        unsigned int key = Utilities::hash_char( "test" );
        if ( key == 2087956275 )
            ut.passes( "Got the expected hash key" );
        else
            ut.failure( "Got the expected hash key" );

        // Test the factor function
        std::vector<int> factors = Utilities::factor( 13958 );
        if ( factors.size() == 3 && factors[0] == 2 && factors[1] == 7 && factors[2] == 997 )
            ut.passes( "Correctly factored 13958" );
        else
            ut.failure( "Correctly factored 13958" );

        // Test getSystemMemory
        size_t system_bytes = Utilities::getSystemMemory();
        std::cout << "Total system bytes = " << system_bytes << std::endl;
        if ( system_bytes > 0 )
            ut.passes( "getSystemMemory" );
        else
            ut.failure( "getSystemMemory" );

        // Test the memory usage
        double t0       = Utilities::time();
        size_t n_bytes1 = Utilities::getMemoryUsage();
        double time1    = Utilities::time() - t0;
        uint64_t *tmp   = new uint64_t[0x100000];
        memset( tmp, 0xAA, 0x100000 * sizeof( uint64_t ) );
        Utilities::nullUse( tmp );
        t0              = Utilities::time();
        size_t n_bytes2 = Utilities::getMemoryUsage();
        double time2    = Utilities::time() - t0;
        delete[] tmp;
        t0              = Utilities::time();
        size_t n_bytes3 = Utilities::getMemoryUsage();
        double time3    = Utilities::time() - t0;
        if ( globalComm.getRank() == 0 ) {
            std::cout << "Number of bytes used for a basic test: " << n_bytes1 << ", " << n_bytes2
                      << ", " << n_bytes3 << std::endl;
            std::cout << "   Time to query: " << time1 * 1e6 << " us, " << time2 * 1e6 << " us, "
                      << time3 * 1e6 << " us" << std::endl;
        }
        if ( n_bytes1 == 0 ) {
            ut.failure( "getMemoryUsage returns 0" );
        } else {
            ut.passes( "getMemoryUsage returns non-zero" );
            if ( n_bytes2 > n_bytes1 ) {
                ut.passes( "getMemoryUsage increases size" );
            } else {
#if defined( USE_MAC )
                ut.expected_failure( "getMemoryUsage does not increase size" );
#else
                ut.failure( "getMemoryUsage increases size" );
#endif
            }
            if ( n_bytes1 == n_bytes3 ) {
                ut.passes( "getMemoryUsage decreases size properly" );
            } else {
#if defined( USE_MAC ) || defined( USE_WINDOWS )
                ut.expected_failure( "getMemoryUsage does not decrease size properly" );
#else
                ut.failure( "getMemoryUsage does not decrease size properly" );
#endif
            }
        }

        // Run large memory test of getMemoryUsage
        if ( system_bytes >= 4e9 && globalComm.getRank() == 0 ) {
            // Test getting the memory usage for 2-4 GB bytes
            // Note: we only run this test on machines with more than 4 GB of memory
            n_bytes1       = Utilities::getMemoryUsage();
            uint64_t *tmp2 = new uint64_t[0x10000001]; // Allocate 2^31+8 bytes
            memset( tmp2, 0xAA, 0x10000001 * sizeof( uint64_t ) );
            Utilities::nullUse( tmp );
            n_bytes2 = Utilities::getMemoryUsage();
            for ( int i = 0; i < 10; i++ ) {
                if ( ( tmp2[rand() % 0x1000000] & 0xFF ) != 0xAA )
                    ut.failure( "Internal error" );
            }
            delete[] tmp2;
            tmp2            = nullptr;
            size_t n_bytes3 = Utilities::getMemoryUsage();
            if ( n_bytes2 > 0x80000000 && n_bytes2 < n_bytes1 + 0x81000000 &&
                 abs_diff( n_bytes1, n_bytes3 ) < 50e3 ) {
                ut.passes( "getMemoryUsage correctly handles 2^31 - 2^32 bytes" );
            } else {
                std::cout << "Memtest 2-4 GB failes: " << n_bytes1 << " " << n_bytes2 << " "
                          << n_bytes3 << std::endl;
                ut.failure( "getMemoryUsage correctly handles 2^31 - 2^32 bytes" );
            }
        }
        if ( system_bytes >= 8e9 && globalComm.getRank() == 0 ) {
            // Test getting the memory usage for > 4 GB bytes
            // Note: we only run this test on machines with more than 8 GB of memory
            n_bytes1       = Utilities::getMemoryUsage();
            size_t size    = 0x20000000;
            uint64_t *tmp2 = new uint64_t[size]; // Allocate 2^31+8 bytes
            if ( tmp == nullptr ) {
                ut.expected_failure( "Unable to allocate variable of size 4 GB" );
            } else {
                memset( tmp2, 0xAA, size * sizeof( uint64_t ) );
                Utilities::nullUse( tmp );
                n_bytes2 = Utilities::getMemoryUsage();
                for ( int i = 0; i < 10; i++ ) {
                    if ( ( tmp2[rand() % size] & 0xFF ) != 0xAA )
                        ut.failure( "Internal error" );
                }
                delete[] tmp2;
                tmp2 = nullptr;
                NULL_USE( tmp2 );
                n_bytes3 = Utilities::getMemoryUsage();
                if ( n_bytes2 > 0x100000000 && n_bytes2 < n_bytes1 + 0x110000000 &&
                     abs_diff( n_bytes1, n_bytes3 ) < 50e3 ) {
                    ut.passes( "getMemoryUsage correctly handles memory > 2^32 bytes" );
                } else {
                    std::cout << "Memtest >4 GB failes: " << n_bytes1 << " " << n_bytes2 << " "
                              << n_bytes3 << std::endl;
                    ut.expected_failure( "getMemoryUsage does not handle memory > 2^32 bytes" );
                }
            }
        }

        // Test getting the current call stack
        double ts1      = Utilities::time();
        auto call_stack = get_call_stack();
        double ts2      = Utilities::time();
        if ( globalComm.getRank() == 0 ) {
            std::cout << "Call stack:" << std::endl;
            for ( auto &elem : call_stack )
                std::cout << "   " << elem.print() << std::endl;
            std::cout << "Time to get call stack: " << ts2 - ts1 << std::endl;
        }
        if ( !call_stack.empty() ) {
            ut.passes( "non empty call stack" );
            bool pass = false;
            if ( call_stack.size() > 1 ) {
                if ( call_stack[1].print().find( "get_call_stack" ) != std::string::npos )
                    pass = true;
            }
            if ( pass )
                ut.passes( "call stack decoded function symbols" );
            else
                ut.expected_failure( "call stack was unable to decode function symbols" );
        } else {
            ut.failure( "non empty call stack" );
        }
        ts1        = Utilities::time();
        auto trace = StackTrace::backtrace();
        ts2        = Utilities::time();
        std::cout << "Time to get backtrace: " << ts2 - ts1 << std::endl;

        // Test getting the symbols
        std::vector<void *> address;
        std::vector<char> type;
        std::vector<std::string> obj;
        int rtn = StackTrace::getSymbols( address, type, obj );
        if ( rtn == 0 && !address.empty() )
            ut.passes( "Read symbols from executable" );

        // Test getting the executable
        std::string exe = StackTrace::getExecutable();
        if ( globalComm.getRank() == 0 )
            std::cout << "Executable: " << exe << std::endl;
        if ( exe.find( "test_Utilities" ) != std::string::npos )
            ut.passes( "getExecutable" );
        else
            ut.failure( "getExecutable" );

        // Test deleting and checking if a file exists
        if ( globalComm.getRank() == 0 ) {
            FILE *fid = fopen( "testDeleteFile.txt", "w" );
            fputs( "Temporary test", fid );
            fclose( fid );
            if ( Utilities::fileExists( "testDeleteFile.txt" ) )
                ut.passes( "File exists" );
            else
                ut.failure( "File exists" );
            Utilities::deleteFile( "testDeleteFile.txt" );
            if ( !Utilities::fileExists( "testDeleteFile.txt" ) )
                ut.passes( "File deleted" );
            else
                ut.failure( "File deleted" );
        }

        // Test creating an empty directory
        Utilities::recursiveMkdir( "." );

        // Finished testing, report the results
        ut.report();
        num_failed = ut.NumFailGlobal();
    }

    // Shutdown
    AMPManager::shutdown();

    // Finished successfully
    return num_failed;
}
