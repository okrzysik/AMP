#include "ampmesh/structured/BoxMeshHelpers.h"
#include "utils/Utilities.h"

#include "math.h"

namespace AMP {
namespace Mesh {
namespace BoxMeshHelpers {


/****************************************************************
* Map the logical coordinates to a circle                       *
****************************************************************/
static inline std::pair<double,double> map_c2p( int method, double xc, double yc )
{
    // map xc > 0 and |yc| < xc ≡ d to (xp,yp) in r=1 using the mapping by:
    //    Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
    //       and Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains", 
    //       SIAM REVIEW, Vol. 50, No. 4, pp.723–752 (2008)
    if ( xc<1e-12 && yc<1e-12 )
        return std::make_pair( 0.0, 0.0 );
    double xp = 0;
    double yp = 0;
    const double sqrt2 = 1.414213562373095;
    if ( method == 1 ) {
        yp = yc/sqrt2;
        xp = sqrt( xc*xc - yp*yp );
    } else if ( method == 2 ) {
        double center = xc/sqrt2 - sqrt( 1.0 - xc*xc/2 );
        yp = yc/sqrt2;
        xp = center + sqrt( 1.0 - yp*yp );
    } else if ( method == 3 ) {
        double D = xc * ( 2 - xc ) / sqrt2;
        double center = D - sqrt( 1.0 - D * D );
        yp = ( 2 - xc ) / sqrt2 * yc;
        xp = center + sqrt( 1.0 - yp*yp );
    } else {
        AMP_ERROR( "Invalid method" );
    }
    return std::make_pair( xp, yp );
}
static inline double root( double a, double b, double c, int r )
{
    double x = 0;
    if ( r == 1 )
        x = (-b+sqrt(b*b-4*a*c))/(2*a);
    else
        x = (-b-sqrt(b*b-4*a*c))/(2*a);
    return x;
}
static inline std::pair<double,double> map_p2c( int method, double xp, double yp )
{
    // Perform the inverse mapping as map_c2p
    if ( xp<1e-12 && yp<1e-12 )
        return std::make_pair( 0.0, 0.0 );
    double xc = 0;
    double yc = 0;
    const double sqrt2 = 1.414213562373095;
    if ( method == 1 ) {
        yc = yp*sqrt2;
        xc = sqrt( xp*xp + yp*yp );
    } else if ( method == 2 ) {
        yc = yp*sqrt2;
        double center = xp - sqrt( 1 - yp*yp );
        xc = root( 2, -2*center*sqrt2, 2*(center*center-1), 1 );
    } else if ( method == 3 ) {
        double center = xp - sqrt( 1 - yp*yp );
        double D = root( 2, -2*center, center*center-1, 1 );
        xc = root( 1, -2, D*sqrt2, 2 );
        yc = yp*sqrt2/(2-xc);
    } else {
        AMP_ERROR( "Invalid method" );
    }
    return std::make_pair( xc, yc );
}
std::pair<double,double> map_logical_circle( double r, int method, double x, double y )
{
    // This maps from a a logically rectangular 3D mesh to a sphere mesh using the mapping by:
    //    Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
    //       and Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains", 
    //       SIAM REVIEW, Vol. 50, No. 4, pp.723–752 (2008)
    const double xc = 2 * x - 1; // Change domain to [-1,1]
    const double yc = 2 * y - 1; // Change domain to [-1,1]
    double xp, yp;
    if ( fabs(xc) >= fabs(yc) ) {
        auto p = map_c2p( method, fabs(xc), fabs(yc) );
        xp = p.first;
        yp = p.second;
    } else {
        auto p = map_c2p( method, fabs(yc), fabs(xc) );
        xp = p.second;
        yp = p.first;
    }
    if ( xc < 0.0 )
        xp = -xp;
    if ( yc < 0.0 )
        yp = -yp;
    return std::make_pair( r*xp, r*yp );
}
std::pair<double,double> map_circle_logical( double r, int method, double x, double y )
{
    // Get the points in the unit circle
    double xp = x/r;
    double yp = y/r;
    // Perform the inverse mapping to [-1,1]
    double xc = 0;
    double yc = 0;
    if ( fabs(xp) >= fabs(yp) ) {
        auto p = map_p2c( method, fabs(xp), fabs(yp) );
        xc = p.first;
        yc = p.second;
    } else {
        auto p = map_p2c( method, fabs(yp), fabs(xp) );
        xc = p.second;
        yc = p.first;
    }
    if ( xp < 0.0 )
        xc = -xc;
    if ( yp < 0.0 )
        yc = -yc;
    // Change domain to [0,1]
    return std::make_pair( 0.5*(xc+1), 0.5*(yc+1) );
}


/****************************************************************
* Helper function to map x,y,z logical coordinates in [0,1]     *
* to x,y,z coordinates in a sphere with radius r                *
****************************************************************/
std::tuple<double,double,double> map_logical_sphere( double r, double x, double y, double z )
{
    // This maps from a a logically rectangular 3D mesh to a sphere mesh using the mapping by:
    //    Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
    //       and Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains", 
    //       SIAM REVIEW, Vol. 50, No. 4, pp.723–752 (2008)
    const double sqrt3 = 1.732050807568877;
    double xc = 2 * x - 1; // Change domain to [-1,1]
    double yc = 2 * y - 1; // Change domain to [-1,1]
    double zc = 2 * z - 1; // Change domain to [-1,1]
    double d  = std::max( std::max( fabs( xc ), fabs( yc ) ), fabs( zc ) );
    double r2 = sqrt( xc * xc + yc * yc + zc * zc );
    r2        = std::max( r2, 1e-10 );
    double x2 = r * d * xc / r2;
    double y2 = r * d * yc / r2;
    double z2 = r * d * zc / r2;
    double w  = d * d;
    x2        = w * x2 + r * ( 1 - w ) * xc / sqrt3;
    y2        = w * y2 + r * ( 1 - w ) * yc / sqrt3;
    z2        = w * z2 + r * ( 1 - w ) * zc / sqrt3;
    return std::make_tuple( x2, y2, z2 );
}


/****************************************************************
* Helper function to map x,y logical coordinates in [0,1]       *
* to x,y,z coordinates on the surface of a sphere               *
****************************************************************/
std::tuple<double,double,double> map_logical_sphere_surface( double R, double x, double y )
{
    // This maps from a a logically rectangular 3D mesh to the surface of a sphere using the mapping by:
    // Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids and Finite
    // GeomType::Volume
    //    Methods for PDEs in Circular and Spherical Domains", SIAM REVIEW, Vol. 50, No. 4, pp.
    //    723–752 (2008)
    double x2 = 2 * x - 1;  // Change domain to [-1,1]
    double x3 = fabs(x2);   // We need to make x go from 1:0:1
    // Map x,y to the unit circle
    auto point = map_logical_circle( 1.0, 3, x3, y );
    double xp = point.first;
    double yp = point.second;
    double zp = sqrt( fabs( 1.0 - ( xp * xp + yp * yp ) ) );
    if ( x2 < 0 )
        zp    = -zp;            // negate z in lower hemisphere
    xp *= R;
    yp *= R;
    zp *= R;
    return std::make_tuple( xp, yp, zp );
}
std::pair<double,double> map_sphere_surface_logical( double R, double x, double y, double z )
{
    // Map from the unit circle to logical
    auto point = map_circle_logical( R, 3, x, y );
    // Change logical coordinates to [0,1]
    double xc = point.first;
    double yc = point.second;
    if ( z < 0 )
        xc = -xc;
    xc = 0.5*(xc+1);
    return std::make_pair( xc, yc );
}


/****************************************************************
* Helper function to map x,y,z logical coordinates in [0,1]     *
* to x,y,z coordinates in a shell with r1 <= r <= r2            *
****************************************************************/
std::tuple<double,double,double> map_logical_shell( double r1, double r2, double x, double y, double z )
{
    // This maps from a a logically rectangular 3D mesh to a shell mesh using the mapping by:
    // Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids and Finite
    // GeomType::Volume
    //    Methods for PDEs in Circular and Spherical Domains", SIAM REVIEW, Vol. 50, No. 4, pp.
    //    723–752 (2008)
    double dr = r2 - r1;
    double Rz = r1 + z * dr;    // radius based on z[0,1]
    return map_logical_sphere_surface( Rz, x, y );
}


} // BoxMeshHelpers namespace
} // Mesh namespace
} // AMP namespace

