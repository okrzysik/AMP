#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

/****************************************************************
 * Map the logical coordinates to a circle                       *
 ****************************************************************/
static inline std::array<double, 2> map_c2p( double xc, double yc )
{
    if ( fabs( xc ) < 1e-12 && fabs( yc ) < 1e-12 )
        return { 0.0, 0.0 };
    if ( fabs( yc ) > fabs( xc ) ) {
        auto [yp, xp] = map_c2p( yc, xc );
        return { xp, yp };
    }
    double scale = std::max( { 1.0, xc, yc } );
    if ( scale > 1.0 ) {
        xc /= scale;
        yc /= scale;
    }
    double xp             = 0;
    double yp             = 0;
    const double invsqrt2 = 0.7071067811865475244;
    double D              = invsqrt2 * xc * ( 2 - xc );
    double center         = D - sqrt( 1.0 - D * D );
    yp                    = invsqrt2 * ( 2 - xc ) * yc;
    xp                    = center + sqrt( 1.0 - yp * yp );
    return { scale * xp, scale * yp };
}
static inline std::array<double, 2> map_p2c( double xp, double yp )
{
    // Perform the inverse mapping as map_c2p
    if ( fabs( xp ) < 1e-12 && fabs( yp ) < 1e-12 )
        return { 0.0, 0.0 };
    if ( fabs( yp ) > fabs( xp ) ) {
        auto [yc, xc] = map_p2c( yp, xp );
        return { xc, yc };
    }
    double scale = std::max( sqrt( xp * xp + yp * yp ), 1.0 );
    if ( scale > 1.0 ) {
        xp /= scale;
        yp /= scale;
    }
    double xc          = 0;
    double yc          = 0;
    const double sqrt2 = 1.4142135623730950488;
    auto z             = xp - sqrt( 1 - yp * yp );
    auto D             = 0.5 * ( z + sqrt( 2 - z * z ) );
    xc                 = 1.0 - sqrt( std::max( 1 - D * sqrt2, 0.0 ) );
    yc                 = yp * sqrt2 / ( 2 - xc );
    return { scale * xc, scale * yc };
}
std::array<double, 2> map_logical_circle( double r, double x, double y )
{
    const double xc = 2 * x - 1; // Change domain to [-1,1]
    const double yc = 2 * y - 1; // Change domain to [-1,1]
    auto [xp, yp]   = map_c2p( fabs( xc ), fabs( yc ) );
    if ( xc < 0.0 )
        xp = -xp;
    if ( yc < 0.0 )
        yp = -yp;
    return { r * xp, r * yp };
}
std::array<double, 2> map_circle_logical( double r, double x, double y )
{
    // Get the points in the unit circle
    double xp = x / r;
    double yp = y / r;
    // Perform the inverse mapping to [-1,1]
    auto [xc, yc] = map_p2c( fabs( xp ), fabs( yp ) );
    if ( xp < 0.0 )
        xc = -xc;
    if ( yp < 0.0 )
        yc = -yc;
    // Change domain to [0,1]
    return { 0.5 * ( xc + 1 ), 0.5 * ( yc + 1 ) };
}


// Test the mapping to/from a logical circle
double test_map_logical_circle( int N, double tol )
{
    auto distance = []( double x, double y, std::array<double, 2> xy ) {
        return sqrt( ( x - xy[0] ) * ( x - xy[0] ) + ( y - xy[1] ) * ( y - xy[1] ) );
    };
    std::random_device rd;
    std::mt19937 gen( rd() );
    bool pass      = true;
    const double r = 2.0;
    std::uniform_real_distribution<> dis( 0, 1 );
    double max_error = 0;
    for ( int i = 0; i < N; i++ ) {
        double x  = dis( gen );
        double y  = dis( gen );
        auto p    = map_logical_circle( r, x, y );
        auto p2   = map_circle_logical( r, p[0], p[1] );
        double r2 = sqrt( p[0] * p[0] + p[1] * p[1] );
        pass      = pass && r2 < r + 1e-15;
        pass      = pass && distance( x, y, p2 ) < tol;
        if ( !( distance( x, y, p2 ) < tol ) )
            printf( "1: %e %e %e %e %e %e\n", x, y, p2[0], p2[1], p2[0] - x, p2[1] - y );
        max_error = std::max( max_error, fabs( distance( x, y, p2 ) ) );
    }
    dis = std::uniform_real_distribution<>( -3.0, 3.0 );
    for ( int i = 0; i < N; i++ ) {
        double x = dis( gen );
        double y = dis( gen );
        auto p   = map_circle_logical( 1.0, x, y );
        auto p2  = map_logical_circle( 1.0, p[0], p[1] );
        pass     = pass && distance( x, y, p2 ) < tol;
        if ( !( distance( x, y, p2 ) < tol ) )
            printf( "2: %e %e %e %e %e %e\n", x, y, p2[0], p2[1], p2[0] - x, p2[1] - y );
        max_error = std::max( max_error, fabs( distance( x, y, p2 ) ) );
    }
    return max_error;
}


// Main function
int main( int, char ** )
{
    double tol = 1e-10;
    double err = test_map_logical_circle( 10000, tol );
    printf( "max error = %e\n", err );
    if ( err < tol ) {
        printf( "Tests passed\n" );
        return 0;
    } else {
        printf( "Tests failed\n" );
        return 1;
    }
}
