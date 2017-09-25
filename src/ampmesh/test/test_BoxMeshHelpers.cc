#include "ampmesh/structured/BoxMeshHelpers.h"

#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#include <math.h>
#include <random>


// Test the mapping to/from a logical circle
void test_map_logical_circle( AMP::UnitTest &ut )
{
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_real_distribution<> dis( 0, 1 );
    const double r = 2.0;
    for ( int method = 1; method <= 3; method++ ) {
        bool pass = true;
        for ( int i = 0; i < 10000; i++ ) {
            double x = dis( gen );
            double y = dis( gen );
            auto p   = AMP::Mesh::BoxMeshHelpers::map_logical_circle( r, method, x, y );
            auto p2 = AMP::Mesh::BoxMeshHelpers::map_circle_logical( r, method, p.first, p.second );
            double r2 = sqrt( p.first * p.first + p.second * p.second );
            pass      = pass && r2 < r + 1e-15;
            pass      = pass && fabs( p2.first - x ) < 1e-10 && fabs( p2.second - y ) < 1e-10;
            if ( fabs( p2.first - x ) > 1e-10 || fabs( p2.second - y ) > 1e-10 )
                printf(
                    "%e %e %e %e %e %e\n", x, y, p2.first, p2.second, p2.first - x, p2.second - y );
        }
        if ( pass )
            ut.passes( "map_logical_circle - " + std::to_string( method ) );
        else
            ut.failure( "map_logical_circle - " + std::to_string( method ) );
    }
}


// Test the mapping to/from the surface of a sphere
void test_map_logical_sphere_surface( AMP::UnitTest &ut )
{
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_real_distribution<> dis( 0, 1 );
    const double r = 2.0;
    bool pass      = true;
    for ( int i = 0; i < 10000; i++ ) {
        double x = dis( gen );
        double y = dis( gen );
        auto p   = AMP::Mesh::BoxMeshHelpers::map_logical_sphere_surface( r, x, y );
        auto p2  = AMP::Mesh::BoxMeshHelpers::map_sphere_surface_logical(
            r, std::get<0>( p ), std::get<1>( p ), std::get<2>( p ) );
        double r2 =
            sqrt( std::get<0>( p ) * std::get<0>( p ) + std::get<1>( p ) * std::get<1>( p ) +
                  std::get<2>( p ) * std::get<2>( p ) );
        pass = pass && r2 < r + 1e-15;
        pass = pass && fabs( p2.first - x ) < 1e-10 && fabs( p2.second - y ) < 1e-10;
        if ( fabs( p2.first - x ) > 1e-10 || fabs( p2.second - y ) > 1e-10 )
            printf( "%e %e %e %e %e %e\n", x, y, p2.first, p2.second, p2.first - x, p2.second - y );
    }
    if ( pass )
        ut.passes( "map_logical_sphere_surface" );
    else
        ut.failure( "map_logical_sphere_surface" );
}


// Main function
int main( int argc, char **argv )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;

    // Run the tests
    test_map_logical_circle( ut );
    test_map_logical_sphere_surface( ut );

    // Print the results and return
    ut.report();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
