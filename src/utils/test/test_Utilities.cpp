#include "UtilityHelpers.h"

#include "AMP/IO/FileSystem.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"

#include "StackTrace/StackTrace.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <vector>


// Run the tests
int main( int argc, char *argv[] )
{

    // Control the behavior of the startup
    AMP::AMPManagerProperties startup_properties;

    // Start AMP
    AMP::AMPManager::startup( argc, argv, startup_properties );
    int num_failed = 0;

    // Limit the scope of variables
    {
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

        // Create the unit test
        AMP::UnitTest ut;

        // Print the banner
        AMP::Utilities::printBanner();

        // Test enable_shared_from_this
        test_shared_from_this( ut );

        // Check the OS
        constexpr auto OS = AMP::Utilities::getOS();
        if constexpr ( OS == AMP::Utilities::OS::Linux )
            ut.passes( "OS: Linux" );
        else if constexpr ( OS == AMP::Utilities::OS::Windows )
            ut.passes( "OS: Windows" );
        else if constexpr ( OS == AMP::Utilities::OS::macOS )
            ut.passes( "OS: macOS" );
        else
            ut.failure( "Unknown OS" );

        // Try converting an int to a string
        if ( AMP::Utilities::intToString( 37, 0 ) == "37" &&
             AMP::Utilities::intToString( 37, 3 ) == "037" )
            ut.passes( "Convert int to string" );
        else
            ut.failure( "Convert int to string" );

        // Test approx_equal
        testApproxEqualInt<int>( ut );
        testApproxEqualInt<unsigned int>( ut );
        testApproxEqualInt<size_t>( ut );
        testApproxEqual<float>( ut );
        testApproxEqual<double>( ut );

        // Test interpolations
        test_interp( ut );

        // Test quicksort performance
        size_t N = 500000;
        std::vector<int> data( N, 31 );
        test_quicksort( ut, data, "identical" );
        for ( size_t i = 0; i < N; i++ )
            data[i] = i;
        test_quicksort( ut, data, "sorted" );
        static std::random_device rd;
        static std::mt19937 gen( rd() );
        static std::uniform_int_distribution<int> dist( 1, 10000000 );
        for ( size_t i = 0; i < N; i++ )
            data[i] = dist( gen );
        test_quicksort( ut, data, "random" );

        // Test quickselect
        for ( size_t i = 0; i < N; i++ )
            data[i] = dist( gen );
        auto data2 = data;
        std::sort( data2.begin(), data2.end() );
        double t    = 0;
        bool pass   = true;
        size_t N_it = 200;
        std::vector<int> data3( data.size(), 0 );
        for ( size_t i = 0; i < N_it; i++ ) {
            data3    = data;
            size_t k = dist( gen ) % N;
            auto t1  = AMP::Utilities::time();
            auto v   = AMP::Utilities::quickselect( data3.size(), data3.data(), k );
            t += AMP::Utilities::time() - t1;
            pass = v == data2[k];
        }
        if ( pass )
            ut.passes( "quickselect" );
        else
            ut.failure( "quickselect" );
        std::cout << "quickselect = " << t / N_it << std::endl;

        // Test the hash key
        unsigned int key = AMP::Utilities::hash_char( "test" );
        if ( key == 2087956275 )
            ut.passes( "Got the expected hash key" );
        else
            ut.failure( "Got the expected hash key" );


        // Test the factor function
        auto factors = AMP::Utilities::factor( 13958 );
        if ( factors == std::vector<int>( { 2, 7, 997 } ) )
            ut.passes( "Correctly factored 13958" );
        else
            ut.failure( "Correctly factored 13958" );
        auto t1 = AMP::Utilities::time();
        N_it    = 10000;
        for ( size_t i = 0; i < N_it; i++ ) {
            [[maybe_unused]] auto tmp = AMP::Utilities::factor( dist( gen ) );
        }
        auto t2 = AMP::Utilities::time();
        std::cout << "factor = " << round( 1e9 * ( t2 - t1 ) / N_it ) << " ns" << std::endl;


        // Test the isPrime function
        if ( !AMP::Utilities::isPrime( 13958 ) && AMP::Utilities::isPrime( 9999991 ) )
            ut.passes( "isPrime" );
        else
            ut.failure( "isPrime" );
        t1 = AMP::Utilities::time();
        for ( size_t i = 0; i < N_it; i++ ) {
            [[maybe_unused]] auto tmp = AMP::Utilities::factor( dist( gen ) );
        }
        t2 = AMP::Utilities::time();
        std::cout << "isPrime = " << round( 1e9 * ( t2 - t1 ) / N_it ) << " ns" << std::endl;


        // Test the primes function
        auto p1 = AMP::Utilities::primes( 50 );
        pass    = p1 ==
               std::vector<uint64_t>( { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47 } );
        t1 = AMP::Utilities::time();
        for ( int i = 0; i < 10; i++ )
            p1 = AMP::Utilities::primes( 1000000 );
        t2 = AMP::Utilities::time();
        if ( p1.size() == 78498 )
            ut.passes( "primes" );
        else
            ut.failure( "primes" );
        std::cout << "size: primes(1000000) = " << p1.size() << std::endl;
        std::cout << "time: primes(1000000) = " << round( 1e6 * ( t2 - t1 ) / 10 ) << " us"
                  << std::endl;


        // Test getting the current call stack
        double ts1      = AMP::Utilities::time();
        auto call_stack = get_call_stack();
        double ts2      = AMP::Utilities::time();
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
        ts1        = AMP::Utilities::time();
        auto trace = StackTrace::backtrace();
        ts2        = AMP::Utilities::time();
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
            if ( AMP::IO::fileExists( "testDeleteFile.txt" ) )
                ut.passes( "File exists" );
            else
                ut.failure( "File exists" );
            AMP::IO::deleteFile( "testDeleteFile.txt" );
            if ( !AMP::IO::fileExists( "testDeleteFile.txt" ) )
                ut.passes( "File deleted" );
            else
                ut.failure( "File deleted" );
        }

        // Test creating directories
        AMP::IO::recursiveMkdir( "." );
        AMP::IO::recursiveMkdir( "testUtilitiesDir/a/b" );
        globalComm.barrier();
        pass = AMP::IO::fileExists( "testUtilitiesDir/a/b" );
        globalComm.barrier();
        if ( globalComm.getRank() == 0 ) {
            AMP::IO::deleteFile( "testUtilitiesDir/a/b" );
            AMP::IO::deleteFile( "testUtilitiesDir/a" );
            AMP::IO::deleteFile( "testUtilitiesDir" );
        }
        globalComm.barrier();
        pass = pass && !AMP::IO::fileExists( "testUtilitiesDir/a/b" );
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

        // Test precision of different types
        test_precision<int32_t>( ut );
        test_precision<int64_t>( ut );
        test_precision<float>( ut );
        test_precision<double>( ut );
        test_precision<long double>( ut );

        // Test printing a warning
        AMP_WARNING( "Testing warning" );

        // Finished testing, report the results
        ut.report();
        num_failed = ut.NumFailGlobal();
    }

    // Shutdown
    AMP::AMPManager::shutdown();

    // Finished successfully
    return num_failed;
}
