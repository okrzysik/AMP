#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

#include "AMP/IO/FileSystem.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/enable_shared_from_this.h"
#include <memory>

#include "StackTrace/StackTrace.h"


using namespace AMP;


// Subract two size_t numbers, returning the absolute value
size_t abs_diff( size_t a, size_t b ) { return ( a >= b ) ? a - b : b - a; }


// This checks approx_equal
template<class T>
void testApproxEqualInt( UnitTest *ut )
{
    std::string type = typeid( T ).name();
    if ( Utilities::approx_equal<T>( 100000, 100000 ) &&
         Utilities::approx_equal_abs<T>( 100000, 100000 ) &&
         !Utilities::approx_equal<T>( 100000, 100001 ) &&
         !Utilities::approx_equal_abs<T>( 100000, 100001 ) )
        ut->passes( "Integer (" + type + ") passes simple check." );
    else
        ut->failure( "Integer (" + type + ") passes simple check." );

    if ( Utilities::approx_equal_abs<T>( 100001, 100000, 1 ) &&
         !Utilities::approx_equal_abs<T>( 100002, 100000, 1 ) )
        ut->passes( "Integer (" + type + ") passes close simple check." );
    else
        ut->failure( "Integer (" + type + ") passes close simple check." );
}
template<class T>
void testApproxEqual( UnitTest *ut )
{
    std::string type = typeid( T ).name();

    T mine      = 1.0;
    T close_rel = mine * static_cast<T>( 1.0 + pow( std::numeric_limits<T>::epsilon(), (T) 0.8 ) );
    T wrong_rel = mine * static_cast<T>( 1.0 + pow( std::numeric_limits<T>::epsilon(), (T) 0.7 ) );
    T close_abs = mine + pow( std::numeric_limits<T>::epsilon(), (T) 0.8 );
    T wrong_abs = mine + pow( std::numeric_limits<T>::epsilon(), (T) 0.7 );
    if ( Utilities::approx_equal( mine, close_rel ) &&
         Utilities::approx_equal_abs( mine, close_abs ) &&
         !Utilities::approx_equal( mine, wrong_rel ) &&
         !Utilities::approx_equal_abs( mine, wrong_abs ) )
        ut->passes( type + " passes simple check near 1" );
    else
        ut->failure( type + " passes simple check near 1" );

    mine      = static_cast<T>( 1e-6 );
    close_rel = mine * static_cast<T>( 1.0 + pow( std::numeric_limits<T>::epsilon(), (T) 0.8 ) );
    wrong_rel = mine * static_cast<T>( 1.0 + pow( std::numeric_limits<T>::epsilon(), (T) 0.7 ) );
    close_abs = mine + pow( std::numeric_limits<T>::epsilon(), (T) 0.8 );
    wrong_abs = mine + pow( std::numeric_limits<T>::epsilon(), (T) 0.7 );
    if ( Utilities::approx_equal( mine, close_rel ) &&
         Utilities::approx_equal_abs( mine, close_abs ) &&
         !Utilities::approx_equal( mine, wrong_rel ) &&
         !Utilities::approx_equal_abs( mine, wrong_abs ) )
        ut->passes( type + " passes simple check near 1e-6" );
    else
        ut->failure( type + " passes simple check near 1e-6" );

    mine      = static_cast<T>( -1e-32 );
    close_rel = mine * static_cast<T>( 1.0 + pow( std::numeric_limits<T>::epsilon(), (T) 0.8 ) );
    wrong_rel = mine * static_cast<T>( 1.0 + pow( std::numeric_limits<T>::epsilon(), (T) 0.7 ) );
    close_abs = mine + pow( std::numeric_limits<T>::epsilon(), (T) 0.8 );
    wrong_abs = mine + pow( std::numeric_limits<T>::epsilon(), (T) 0.7 );
    if ( Utilities::approx_equal( mine, close_rel ) &&
         Utilities::approx_equal_abs( mine, close_abs ) &&
         !Utilities::approx_equal( mine, wrong_rel ) &&
         !Utilities::approx_equal_abs( mine, wrong_abs ) )
        ut->passes( type + " passes simple check near -1e-32" );
    else
        ut->failure( type + " passes simple check near -1e-32" );
}


// Function to return the call stack
std::vector<StackTrace::stack_info> get_call_stack() { return StackTrace::getCallStack(); }


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
    std::shared_ptr<dummy> getPtr() { return shared_from_this(); }
};
static inline bool test_shared_from_this_pointer( const std::shared_ptr<dummy> &p1 )
{
    bool pass = p1.use_count() == 1;
    int count1, count2;
    auto p2 = p1->getPtr();
    count1  = p1.use_count();
    count2  = p2.use_count();
    pass    = pass && count1 == 2 && count2 == 2;
    std::shared_ptr<dummy> p3( p1 );
    count1 = p2.use_count();
    count2 = p3.use_count();
    pass   = pass && count1 == 3 && count2 == 3;
    p2.reset();
    count1  = p2.use_count();
    count2  = p3.use_count();
    pass    = pass && count1 == 0 && count2 == 2;
    auto p4 = p3->getPtr();
    count1  = p3.use_count();
    count2  = p4.use_count();
    pass    = pass && count1 == 3 && count2 == 3;
    std::shared_ptr<dummy> p5( p3.get(), []( void * ) {} );
    count1 = p3.use_count();
    count2 = p5.use_count();
    pass   = pass && count1 == 3 && count2 == 1;
    p5.reset();
    count1 = p3.use_count();
    count2 = p5.use_count();
    pass   = pass && p3.use_count() == 3 && p5.use_count() == 0;
    return pass;
}
void test_shared_from_this( UnitTest *ut )
{
    bool pass = true;
    try {
        auto ptr = std::make_shared<dummy>();
        pass     = test_shared_from_this_pointer( ptr );
    } catch ( ... ) {
        pass = false;
    }
    if ( pass )
        ut->passes( "shared_from_this 1" );
    else
        ut->failure( "shared_from_this 1" );
    try {
        auto *p1 = new dummy;
        auto p2  = p1->getPtr();
        pass     = test_shared_from_this_pointer( p2 );
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

        // Check the OS
        constexpr auto OS = Utilities::getOS();
        if constexpr ( OS == Utilities::OS::Linux )
            ut.passes( "OS: Linux" );
        else if constexpr ( OS == Utilities::OS::Windows )
            ut.passes( "OS: Windows" );
        else if constexpr ( OS == Utilities::OS::macOS )
            ut.passes( "OS: macOS" );
        else
            ut.failure( "Known OS" );

        // Try converting an int to a string
        if ( Utilities::intToString( 37, 0 ) == "37" && Utilities::intToString( 37, 3 ) == "037" )
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
        size_t N = 100000;
        std::vector<int> data1( N );
        srand( static_cast<unsigned int>( time( nullptr ) ) );
        for ( size_t i = 0; i < N; i++ )
            data1[i] = rand();
        auto data2 = data1;
        auto data3 = data1;
        auto data4 = data1;
        double t1  = Utilities::time();
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

        // Test quickselect
        for ( size_t i = 0; i < N; i++ )
            data1[i] = rand();
        data2 = data1;
        std::sort( data2.begin(), data2.end() );
        double t    = 0;
        pass        = true;
        size_t N_it = 200;
        for ( size_t i = 0; i < N_it; i++ ) {
            data3    = data1;
            size_t k = rand() % N;
            t1       = Utilities::time();
            auto v   = Utilities::quickselect( data3.size(), data3.data(), k );
            t += Utilities::time() - t1;
            pass = v == data2[k];
        }
        if ( pass )
            ut.passes( "quickselect" );
        else
            ut.failure( "quickselect" );
        std::cout << "quickselect = " << t / N_it << std::endl;

        // Test the hash key
        unsigned int key = Utilities::hash_char( "test" );
        if ( key == 2087956275 )
            ut.passes( "Got the expected hash key" );
        else
            ut.failure( "Got the expected hash key" );


        // Test the factor function
        auto factors = Utilities::factor( 13958 );
        if ( factors == std::vector<int>( { 2, 7, 997 } ) )
            ut.passes( "Correctly factored 13958" );
        else
            ut.failure( "Correctly factored 13958" );
        std::default_random_engine gen;
        std::uniform_int_distribution<int> dist( 1, 10000000 );
        t1   = Utilities::time();
        N_it = 10000;
        for ( size_t i = 0; i < N_it; i++ ) {
            auto tmp = AMP::Utilities::factor( dist( gen ) );
            NULL_USE( tmp );
        }
        t2 = Utilities::time();
        std::cout << "factor = " << round( 1e9 * ( t2 - t1 ) / N_it ) << " ns" << std::endl;


        // Test the isPrime function
        if ( !Utilities::isPrime( 13958 ) && Utilities::isPrime( 9999991 ) )
            ut.passes( "isPrime" );
        else
            ut.failure( "isPrime" );
        t1 = Utilities::time();
        for ( size_t i = 0; i < N_it; i++ ) {
            auto tmp = AMP::Utilities::factor( dist( gen ) );
            NULL_USE( tmp );
        }
        t2 = Utilities::time();
        std::cout << "isPrime = " << round( 1e9 * ( t2 - t1 ) / N_it ) << " ns" << std::endl;


        // Test the primes function
        auto p1 = Utilities::primes( 50 );
        pass    = p1 ==
               std::vector<uint64_t>( { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47 } );
        t1 = Utilities::time();
        for ( int i = 0; i < 10; i++ )
            p1 = Utilities::primes( 1000000 );
        t2 = Utilities::time();
        if ( p1.size() == 78498 )
            ut.passes( "primes" );
        else
            ut.failure( "primes" );
        std::cout << "size: primes(1000000) = " << p1.size() << std::endl;
        std::cout << "time: primes(1000000) = " << round( 1e6 * ( t2 - t1 ) / 10 ) << " us"
                  << std::endl;


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
            pass = false;
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
            if ( IO::fileExists( "testDeleteFile.txt" ) )
                ut.passes( "File exists" );
            else
                ut.failure( "File exists" );
            IO::deleteFile( "testDeleteFile.txt" );
            if ( !IO::fileExists( "testDeleteFile.txt" ) )
                ut.passes( "File deleted" );
            else
                ut.failure( "File deleted" );
        }

        // Test creating directories
        IO::recursiveMkdir( "." );
        IO::recursiveMkdir( "testUtilitiesDir/a/b" );
        globalComm.barrier();
        pass = IO::fileExists( "testUtilitiesDir/a/b" );
        if ( globalComm.getRank() == 0 ) {
            IO::deleteFile( "testUtilitiesDir/a/b" );
            IO::deleteFile( "testUtilitiesDir/a" );
            IO::deleteFile( "testUtilitiesDir" );
        }
        globalComm.barrier();
        pass = pass && !IO::fileExists( "testUtilitiesDir/a/b" );
        if ( pass )
            ut.passes( "Create/destroy directory" );
        else
            ut.failure( "Create/destroy directory" );

        // Test catching an error
        try {
            AMP_ERROR( "test_error" );
            ut.failure( "Failed to catch error" );
        } catch ( const StackTrace::abort_error &err ) {
            if ( err.message == "test_error" )
                ut.passes( "Caught error" );
            else
                ut.failure( "Failed to catch error with proper message" );
        } catch ( std::exception &err ) {
            ut.failure( "Caught unknown exception type" );
        }

        // Test printing a warning
        AMP_WARNING( "Testing warning" );

        // Finished testing, report the results
        ut.report();
        num_failed = ut.NumFailGlobal();
    }

    // Shutdown
    AMPManager::shutdown();

    // Finished successfully
    return num_failed;
}
