#include "ampmesh/structured/BoxMeshHelpers.h"
#include "utils/Utilities.h"

#include "math.h"

namespace AMP {
namespace Mesh {
namespace BoxMeshHelpers {


/****************************************************************
* Map the logical coordinates to a circle                       *
****************************************************************/
std::pair<double,double> map_logical_circle( double r, int method, double x, double y )
{
    // This maps from a a logically rectangular 3D mesh to a sphere mesh using the mapping by:
    //    Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
    //       and Finite Volume Methods for PDEs in Circular and Spherical Domains", 
    //       SIAM REVIEW, Vol. 50, No. 4, pp.723–752 (2008)
    const double sqrt2 = 1.41421356237;
    // map [0,1] x [0,1] to circle of radius r
    double xc = 2 * x - 1; // Change domain to [-1,1]
    double yc = 2 * y - 1; // Change domain to [-1,1]
    if ( fabs( xc ) < 1e-12 && fabs( yc ) < 1e-12 ) {
        // We are dealing with the center point
        return std::pair<double,double>(0,0);
    }
    double d = std::max( fabs( xc ), fabs( yc ) ); // value on diagonal of computational grid
    double D = 0, R = 0;
    if ( method == 1 ) {
        D = r * d / sqrt2; // mapping d to D(d)
        R = r * d;         // mapping d to R(d)
    } else if ( method == 2 ) {
        D = r * d / sqrt2; // mapping d to D(d)
        R = r;             // mapping d to R(d)
    } else if ( method == 3 ) {
        D = r * d * ( 2 - d ) / sqrt2; // mapping d to D(d)
        R = r;                         // mapping d to R(d)
    } else {
        AMP_ERROR( "Invalid method" );
    }
    double center = D - sqrt( R * R - D * D );
    double xp     = D / d * fabs( xc );
    double yp     = D / d * fabs( yc );
    if ( fabs( yc ) >= fabs( xc ) )
        yp = center + sqrt( R * R - xp * xp );
    if ( fabs( xc ) >= fabs( yc ) )
        xp = center + sqrt( R * R - yp * yp );
    if ( xc < 0.0 )
        xp = -xp;
    if ( yc < 0.0 )
        yp = -yp;
    return std::make_pair( xp, yp );
}


/****************************************************************
* Helper function to map x,y,z logical coordinates in [0,1]     *
* to x,y,z coordinates in a sphere with radius r                *
****************************************************************/
std::tuple<double,double,double> map_logical_sphere( double r, double x, double y, double z )
{
    // This maps from a a logically rectangular 3D mesh to a sphere mesh using the mapping by:
    //    Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
    //       and Finite Volume Methods for PDEs in Circular and Spherical Domains", 
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
* Helper function to map x,y,z logical coordinates in [0,1]     *
* to x,y,z coordinates in a shell with r1 <= r <= r2            *
****************************************************************/
std::tuple<double,double,double> map_logical_shell( double r1, double r2, double x, double y, double z )
{
    // This maps from a a logically rectangular 3D mesh to a shell mesh using the mapping by:
    // Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids and Finite
    // Volume
    //    Methods for PDEs in Circular and Spherical Domains", SIAM REVIEW, Vol. 50, No. 4, pp.
    //    723–752 (2008)
    double dr = r2 - r1;
    double x2 = 2 * x - 1;  // Change domain to [-1,1]
    double x3 = fabs(x2);   // We need to make x go from 1:0:1
    // Map x,y to the unit circle
    auto point = map_logical_circle( 1.0, 3, x3, y );
    double xp = point.first;
    double yp = point.second;
    double zp = sqrt( fabs( 1.0 - ( xp * xp + yp * yp ) ) );
    if ( x2 < 0 )
        zp    = -zp;            // negate z in lower hemisphere
    double Rz = r1 + z * dr;    // radius based on z[0,1]
    xp *= Rz;
    yp *= Rz;
    zp *= Rz;
    return std::make_tuple( xp, yp, zp );
}


} // BoxMeshHelpers namespace
} // Mesh namespace
} // AMP namespace

