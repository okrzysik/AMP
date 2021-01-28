#include "AMP/ampmesh/MeshPoint.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"

#include <cmath>
#include <random>


using namespace AMP::Geometry::GeometryHelpers;
using AMP::Mesh::Point;


// Test the mapping to/from a logical circle
void test_dist_line( AMP::UnitTest &ut )
{
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_real_distribution<> dis( -50, 50 );
    // Test ray - line-segment intersection
    bool pass = true;
    for ( int i = 0; i < 1000; i++ ) {
        // Generate an intersection point, then the rays and line segments
        Point pi  = { dis( gen ), dis( gen ) };
        Point p0  = { dis( gen ), dis( gen ) };
        Point p1  = p0 + 1.34 * ( pi - p0 );
        Point p2  = p0 + 0.98 * ( pi - p0 );
        Point p   = { dis( gen ), dis( gen ) };
        Point dir = normalize( pi - p );
        double d  = ( pi.x() - p.x() ) / dir.x();
        double d1 = distanceToLine( p, dir, p0, p1 );
        double d2 = distanceToLine( p, dir, p0, p2 );
        double d3 = distanceToLine( p, -dir, p0, p2 );
        double d4 = distanceToLine( p, dir, p1, p0 );
        pass      = pass && fabs( d - d1 ) < 1e-10;
        pass      = pass && d2 == std::numeric_limits<double>::infinity();
        pass      = pass && d3 == std::numeric_limits<double>::infinity();
        pass      = pass && fabs( d - d4 ) < 1e-10;
        if ( !( fabs( d - d1 ) < 1e-10 ) )
            printf( "distanceToLine: %f %f\n", d, d1 );
    }
    if ( pass )
        ut.passes( "distanceToLine" );
    else
        ut.failure( "distanceToLine" );
}


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
            double x  = dis( gen );
            double y  = dis( gen );
            auto p    = map_logical_circle( r, method, x, y );
            auto p2   = map_circle_logical( r, method, p.first, p.second );
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


// Test the mapping to/from a regular polygon
void test_map_logical_poly( AMP::UnitTest &ut )
{
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_real_distribution<> dis( 0, 1 );
    const double r = 2.3;
    for ( int N = 3; N <= 10; N++ ) {
        bool pass = true;
        for ( int i = 0; i < 10000; i++ ) {
            double x  = dis( gen );
            double y  = dis( gen );
            auto p    = map_logical_poly( N, r, x, y );
            auto p2   = map_poly_logical( N, r, p.first, p.second );
            double r2 = sqrt( p.first * p.first + p.second * p.second );
            pass      = pass && r2 < r + 1e-15;
            pass      = pass && fabs( p2.first - x ) < 1e-10 && fabs( p2.second - y ) < 1e-10;
            if ( fabs( p2.first - x ) > 1e-10 || fabs( p2.second - y ) > 1e-10 )
                printf(
                    "%e %e %e %e %e %e\n", x, y, p2.first, p2.second, p2.first - x, p2.second - y );
        }
        if ( pass )
            ut.passes( "map_logical_poly - " + std::to_string( N ) );
        else
            ut.failure( "map_logical_poly - " + std::to_string( N ) );
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
        double x  = dis( gen );
        double y  = dis( gen );
        auto p    = map_logical_sphere_surface( r, x, y );
        auto p2   = map_sphere_surface_logical( r, p[0], p[1], p[2] );
        double r2 = sqrt( p[0] * p[0] + p[1] * p[1] + p[2] * p[2] );
        pass      = pass && r2 < r + 1e-15;
        pass      = pass && fabs( p2.first - x ) < 1e-10 && fabs( p2.second - y ) < 1e-10;
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
    test_dist_line( ut );
    test_map_logical_poly( ut );
    test_map_logical_circle( ut );
    test_map_logical_sphere_surface( ut );

    // Print the results and return
    ut.report();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
