#include <iostream>
#include <math.h>
#include <stdexcept>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <vector>

#include "utils/AMPManager.h"
#include "utils/Array.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

// The main function
int main( int argc, char *argv[] )
{
    // Startup
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;

    // these are currently not defined for AMP
    //    AMP::Utilities::setAbortBehavior(true,true,true);
    //    AMP::Utilities::setErrorHandlers();

    // Limit the scope of variables
    {
        // Create several matrices
        AMP::Array<double> M1, M2( 10, 5 );
        M1.resize( 10, 7 );
        for ( size_t i = 0; i < M2.size( 0 ); i++ ) {
            for ( size_t j = 0; j < M2.size( 1 ); j++ ) {
                M1( i, j ) = i + 10 * j;
                M2( i, j ) = i + 10 * j;
            }
        }
        M1.resize( 10, 5 );
        AMP::Array<double> M3( M1 );
        AMP::Array<double> M4 = M2;
        AMP::Array<double> M5 = M1;
        M5( 0, 0 ) = -1;
        if ( M1 == M2 && M1 == M3 && M1 == M4 && M1 != M5 )
            ut.passes( "Array constructors" );
        else
            ut.failure( "Array constructors" );
        // Test std::string
        bool pass = true;
        AMP::Array<std::string> S;
        pass = pass && S.length() == 0;
        S.resize( 1 );
        pass   = pass && S.length() == 1;
        pass   = pass && S( 0 ).size() == 0;
        S( 0 ) = std::string( "test" );
        pass   = pass && S( 0 ) == "test";
        if ( pass )
            ut.passes( "AMP::Array string" );
        else
            ut.failure( "AMP::Array string" );
        // Test a failed allocation
        try {
            size_t N = 10000;
            AMP::Array<double> M( N, N, N );
#if defined( __APPLE__ )
            ut.expected_failure( "Failed allocation succeeded (MAC)" );
#else
            ut.failure( "Failed allocation succeeded???" );
#endif
            AMP_ASSERT( M.length() == N * N * N );
        } catch ( ... ) {
            ut.passes( "Caught failed allocation" );
        }
        // Test math operators
        if ( M1.min() == 0 )
            ut.passes( "min" );
        else
            ut.failure( "min" );
        if ( M1.max() == 49 )
            ut.passes( "max" );
        else
            ut.failure( "max" );
        if ( M1.sum() == 1225 )
            ut.passes( "sum" );
        else
            ut.failure( "sum" );
        if ( M1.mean() == 24.5 )
            ut.passes( "mean" );
        else
            ut.failure( "mean" );
        if ( !M1.NaNs() )
            ut.passes( "NaNs" );
        else
            ut.failure( "NaNs" );
        // Test math operators with index subsets
        std::vector<size_t> idx{ 0, 4, 0, 2 };
        if ( M1.min( idx ) == 0 )
            ut.passes( "min on subset" );
        else
            ut.failure( "min on subset" );
        if ( M1.max( idx ) == 24 )
            ut.passes( "max on subset" );
        else
            ut.failure( "max on subset" );
        if ( M1.sum( idx ) == 180 )
            ut.passes( "sum on subset" );
        else {
            ut.failure( "sum on subset" );
        }
        if ( M1.mean( idx ) == 12 )
            ut.passes( "mean on subset" );
        else
            ut.failure( "mean on subset" );
        // Test find
        std::vector<size_t> index = M1.find( 7, []( double a, double b ) { return a == b; } );
        if ( index.size() != 1 )
            ut.failure( "find" );
        else if ( index[0] == 7 )
            ut.passes( "find" );
        else
            ut.failure( "find" );
        // Test subset
        M3 = M1.subset( { 0, 9, 0, 4 } );
        if ( M3 == M1 )
            ut.passes( "full subset" );
        else
            ut.failure( "full subset" );
        M3   = M1.subset( { 3, 7, 1, 3 } );
        pass = true;
        for ( size_t i = 0; i < M3.size( 0 ); i++ ) {
            for ( size_t j = 0; j < M3.size( 1 ); j++ )
                pass = pass && M3( i, j ) == ( i + 3 ) + 10 * ( j + 1 );
        }
        if ( pass )
            ut.passes( "partial subset" );
        else
            ut.failure( "partial subset" );
        M3.scale( 2 );
        M2.copySubset( { 3, 7, 1, 3 }, M3 );
        pass = true;
        for ( size_t i = 0; i < M3.size( 0 ); i++ ) {
            for ( size_t j = 0; j < M3.size( 1 ); j++ )
                pass = pass && M3( i, j ) == M2( i + 3, j + 1 );
        }
        if ( pass )
            ut.passes( "copyFromSubset" );
        else
            ut.failure( "copyFromSubset" );
        // Test the time required to create a view
        AMP::Array<double> M_view;
        double t1 = AMP::Utilities::time();
        for ( size_t i = 0; i < 100000; i++ ) {
            M_view.viewRaw( { M1.size( 0 ), M1.size( 1 ) }, M1.data() );
            NULL_USE( M_view );
        }
        double t2 = AMP::Utilities::time();
        if ( M_view == M1 )
            ut.passes( "view" );
        else
            ut.failure( "view" );
        AMP::pout << "Time to create view: " << ( t2 - t1 ) * 1e9 / 100000 << " ns\n";
    }

    // Finished
    ut.report();
    int num_failed = static_cast<int>( ut.NumFailGlobal() );
    if ( num_failed == 0 )
        AMP::pout << "All tests passed\n";
    AMP::AMPManager::shutdown();
    return num_failed;
}
