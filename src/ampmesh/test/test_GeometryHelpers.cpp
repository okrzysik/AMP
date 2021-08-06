#include "AMP/ampmesh/MeshPoint.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <random>


using namespace AMP::Geometry::GeometryHelpers;
using AMP::Mesh::Point;


/****************************************************************
 * Vector operations                                             *
 ****************************************************************/
static inline Point3D operator+( const Point3D &x, const Point3D &y )
{
    return { x[0] + y[0], x[1] + y[1], x[2] + y[2] };
}
static inline Point3D operator-( const Point3D &x, const Point3D &y )
{
    return { x[0] - y[0], x[1] - y[1], x[2] - y[2] };
}
static inline Point3D operator*( double x, const Point3D &y )
{
    return { x * y[0], x * y[1], x * y[2] };
}
static inline double dot( const Point3D &x, const Point3D &y )
{
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}
static inline Point3D cross( const Point3D &x, const Point3D &y )
{
    return { x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0] };
}
static inline Point3D normalize( const Point3D &x )
{
    double tmp = 1.0 / sqrt( dot( x, x ) );
    return { tmp * x[0], tmp * x[1], tmp * x[2] };
}


// Test the mapping to/from a logical circle
void test_dist_line( int N, AMP::UnitTest &ut )
{
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_real_distribution<> dis( -50, 50 );
    // Test ray - line-segment intersection
    bool pass = true;
    auto t1   = std::chrono::high_resolution_clock::now();
    for ( int i = 0; i < N; i++ ) {
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
        pass      = pass && fabs( d - d1 ) < 1e-8;
        pass      = pass && d2 == std::numeric_limits<double>::infinity();
        pass      = pass && d3 == std::numeric_limits<double>::infinity();
        pass      = pass && fabs( d - d4 ) < 1e-8;
        if ( !( fabs( d - d1 ) < 1e-8 ) )
            printf( "distanceToLine: %f %f %e\n", d, d1, d - d1 );
    }
    auto t2    = std::chrono::high_resolution_clock::now();
    int64_t ns = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();
    printf( "distanceToLine: %i ns\n", static_cast<int>( ns / ( 4 * N ) ) );
    if ( pass )
        ut.passes( "distanceToLine" );
    else
        ut.failure( "distanceToLine" );
}


// Test the mapping to/from a logical circle
void test_map_logical_circle( int N, AMP::UnitTest &ut )
{
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_real_distribution<> dis( 0, 1 );
    const double r = 2.0;
    for ( int method = 1; method <= 3; method++ ) {
        bool pass = true;
        auto t1   = std::chrono::high_resolution_clock::now();
        for ( int i = 0; i < N; i++ ) {
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
        auto t2    = std::chrono::high_resolution_clock::now();
        int64_t ns = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();
        printf( "map_logical_circle - %i: %i ns\n", method, static_cast<int>( ns / ( 4 * N ) ) );
        if ( pass )
            ut.passes( "map_logical_circle - " + std::to_string( method ) );
        else
            ut.failure( "map_logical_circle - " + std::to_string( method ) );
    }
}


// Test the mapping to/from a regular polygon
void test_map_logical_poly( int N, AMP::UnitTest &ut )
{
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_real_distribution<> dis( 0, 1 );
    const double r = 2.3;
    for ( int Np = 3; Np <= 10; Np++ ) {
        bool pass = true;
        auto t1   = std::chrono::high_resolution_clock::now();
        for ( int i = 0; i < N; i++ ) {
            double x  = dis( gen );
            double y  = dis( gen );
            auto p    = map_logical_poly( Np, r, x, y );
            auto p2   = map_poly_logical( Np, r, p.first, p.second );
            double r2 = sqrt( p.first * p.first + p.second * p.second );
            pass      = pass && r2 < r + 1e-15;
            pass      = pass && fabs( p2.first - x ) < 1e-10 && fabs( p2.second - y ) < 1e-10;
            if ( fabs( p2.first - x ) > 1e-10 || fabs( p2.second - y ) > 1e-10 )
                printf(
                    "%e %e %e %e %e %e\n", x, y, p2.first, p2.second, p2.first - x, p2.second - y );
        }
        auto t2    = std::chrono::high_resolution_clock::now();
        int64_t ns = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();
        printf( "map_logical_poly - %i: %i ns\n", Np, static_cast<int>( ns / N ) );
        if ( pass )
            ut.passes( "map_logical_poly - " + std::to_string( Np ) );
        else
            ut.failure( "map_logical_poly - " + std::to_string( Np ) );
    }
}


// Test the mapping to/from the surface of a sphere
void test_map_logical_sphere_surface( int N, AMP::UnitTest &ut )
{
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_real_distribution<> dis( 0, 1 );
    const double r = 2.0;
    bool pass      = true;
    auto t1        = std::chrono::high_resolution_clock::now();
    for ( int i = 0; i < N; i++ ) {
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
    auto t2    = std::chrono::high_resolution_clock::now();
    int64_t ns = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();
    printf( "map_logical_sphere_surface: %i ns\n", static_cast<int>( ns / N ) );
    if ( pass )
        ut.passes( "map_logical_sphere_surface" );
    else
        ut.failure( "map_logical_sphere_surface" );
}


// Create random triangle/ray/distance sets in the same plane
static std::tuple<std::array<Point3D, 3>, Point3D, Point3D, double> createTriRayPlane()
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_real_distribution<> dis( -1.0, 1.0 );
    // Create the center of a triangle and point of intersection
    Point3D c  = { dis( gen ), dis( gen ), dis( gen ) };
    Point3D pi = { dis( gen ), dis( gen ), dis( gen ) };
    // Create the ray
    auto dir = normalize( c - pi );
    double d = 1.0 + dis( gen );
    auto pos = pi - d * dir;
    // Create the triangle points
    Point3D v1 = normalize( { dis( gen ), dis( gen ), dis( gen ) } );
    std::array<Point3D, 3> tri;
    tri[0] = pi + fabs( dis( gen ) ) * v1;
    tri[1] = pi - fabs( dis( gen ) ) * v1;
    tri[2] = pi + 2 * dir;
    return std::tie( tri, pos, dir, d );
}


// Create random triangle/ray/distance sets in the same plane
static std::tuple<std::array<Point3D, 3>, Point3D, Point3D, double> createTriRayVol()
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_real_distribution<> dis( -1.0, 1.0 );
    // Create the triangle points
    std::array<Point3D, 3> tri;
    tri[0] = { dis( gen ), dis( gen ), dis( gen ) };
    tri[1] = { dis( gen ), dis( gen ), dis( gen ) };
    tri[2] = { dis( gen ), dis( gen ), dis( gen ) };
    // Create a random point within the triangle and direction
    Point3D L = { fabs( dis( gen ) ), fabs( dis( gen ) ), fabs( dis( gen ) ) };
    double t  = 1.0 / ( L[0] + L[1] + L[2] );
    L         = { t * L[0], t * L[1], t * L[2] };
    auto pi   = L[0] * tri[0] + L[1] * tri[1] + L[2] * tri[2];
    auto dir  = normalize( { dis( gen ), dis( gen ), dis( gen ) } );
    double d  = 1.0 + fabs( dis( gen ) );
    auto pos  = pi - d * dir;
    return std::tie( tri, pos, dir, d );
}


// Test ray-triangle intersection
void test_ray_triangle_intersection( int N, AMP::UnitTest &ut )
{
    bool pass = true;
    // Create some sample triangles
    std::array<Point2D, 3> t1 = { { { 0, 0 }, { 1, 0 }, { 0, 1 } } };
    std::array<Point3D, 3> t2 = { { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 } } };
    // Test 2D triangle
    double d11 = distanceToTriangle( t1, { -1, 0 }, { 1, 0 } );
    double d12 = distanceToTriangle( t1, { -1, 0.5 }, { 1, 0 } );
    double d13 = distanceToTriangle( t1, { -1, -0.5 }, { 1, 0 } );
    double d14 = distanceToTriangle( t1, { 0.1, 0.2 }, { -1, 0 } );
    pass       = pass && fabs( d11 - 1.0 ) < 1e-12;
    pass       = pass && fabs( d12 - 1.0 ) < 1e-12;
    pass       = pass && d13 == std::numeric_limits<double>::infinity();
    pass       = pass && fabs( d14 + 0.1 ) < 1e-12;
    // Test 2D triangles in 3D
    double d21 = distanceToTriangle( t2, { -1, 0, 0 }, { 1, 0, 0 } );
    double d22 = distanceToTriangle( t2, { -1, 0.5, 0 }, { 1, 0, 0 } );
    double d23 = distanceToTriangle( t2, { -1, -0.5, 0 }, { 1, 0, 0 } );
    pass       = pass && fabs( d21 - 1.0 ) < 1e-12;
    pass       = pass && fabs( d22 - 1.0 ) < 1e-12;
    pass       = pass && d23 == std::numeric_limits<double>::infinity();
    auto start = std::chrono::high_resolution_clock::now();
    for ( int i = 0; i < N; i++ ) {
        auto [tri, pos, dir, d] = createTriRayPlane();
        double d2               = distanceToTriangle( tri, pos, dir );
        if ( fabs( d - d2 ) > 1e-8 ) {
            pass = false;
        }
    }
    for ( int i = 0; i < N; i++ ) {
        auto [tri, pos, dir, d] = createTriRayVol();
        double d2               = distanceToTriangle( tri, pos, dir );
        if ( fabs( d - d2 ) > 1e-8 ) {
            pass = false;
        }
    }
    auto end   = std::chrono::high_resolution_clock::now();
    int64_t ns = std::chrono::duration_cast<std::chrono::nanoseconds>( end - start ).count();
    printf( "ray-triangle intersection: %i ns\n", static_cast<int>( ns / ( 2 * N ) ) );
    if ( pass ) {
        ut.passes( "ray-triangle intersection" );
    } else {
        ut.failure( "ray-triangle intersection" );
    }
}


// Main function
int main( int argc, char **argv )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;

    // Run the tests
    test_dist_line( 10000, ut );
    test_map_logical_poly( 10000, ut );
    test_map_logical_circle( 10000, ut );
    test_map_logical_sphere_surface( 10000, ut );
    test_ray_triangle_intersection( 1000, ut );

    // Print the results and return
    ut.report();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
