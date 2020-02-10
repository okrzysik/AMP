#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/ampmesh/Geometry.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>
#include <cmath>

namespace AMP {
namespace Geometry {
namespace GeometryHelpers {


/****************************************************************
 * Map the logical coordinates to a circle                       *
 ****************************************************************/
static inline std::pair<double, double> map_c2p( int method, double xc, double yc )
{
    // map xc > 0 and |yc| < xc ≡ d to (xp,yp) in r=1 using the mapping by:
    //    Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
    //       and Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains",
    //       SIAM REVIEW, Vol. 50, No. 4, pp.723–752 (2008)
    if ( xc < 1e-12 && yc < 1e-12 )
        return std::make_pair( 0.0, 0.0 );
    double xp          = 0;
    double yp          = 0;
    const double sqrt2 = 1.414213562373095;
    if ( method == 1 ) {
        yp = yc / sqrt2;
        xp = sqrt( xc * xc - yp * yp );
    } else if ( method == 2 ) {
        double center = xc / sqrt2 - sqrt( 1.0 - xc * xc / 2 );
        yp            = yc / sqrt2;
        xp            = center + sqrt( 1.0 - yp * yp );
    } else if ( method == 3 ) {
        double D      = xc * ( 2 - xc ) / sqrt2;
        double center = D - sqrt( 1.0 - D * D );
        yp            = ( 2 - xc ) / sqrt2 * yc;
        xp            = center + sqrt( 1.0 - yp * yp );
    } else {
        AMP_ERROR( "Invalid method" );
    }
    return std::make_pair( xp, yp );
}
static inline double root( double a, double b, double c, int r )
{
    double t = b * b - 4 * a * c;
    if ( fabs( t ) < 1e-12 * fabs( b * b ) )
        return -b / ( 2 * a );
    double s = r == 1 ? 1 : -1;
    double x = ( -b + s * sqrt( t ) ) / ( 2 * a );
    return x;
}
static inline std::pair<double, double> map_p2c( int method, double xp, double yp )
{
    // Perform the inverse mapping as map_c2p
    if ( xp < 1e-12 && yp < 1e-12 )
        return std::make_pair( 0.0, 0.0 );
    double xc          = 0;
    double yc          = 0;
    const double sqrt2 = 1.414213562373095;
    if ( method == 1 ) {
        yc = yp * sqrt2;
        xc = sqrt( xp * xp + yp * yp );
    } else if ( method == 2 ) {
        yc            = yp * sqrt2;
        double center = xp - sqrt( 1 - yp * yp );
        xc            = root( 2, -2 * center * sqrt2, 2 * ( center * center - 1 ), 1 );
    } else if ( method == 3 ) {
        double center = xp - sqrt( 1 - yp * yp );
        double D      = root( 2, -2 * center, center * center - 1, 1 );
        xc            = root( 1, -2, D * sqrt2, 2 );
        yc            = yp * sqrt2 / ( 2 - xc );
    } else {
        AMP_ERROR( "Invalid method" );
    }
    return std::make_pair( xc, yc );
}
std::pair<double, double> map_logical_circle( double r, int method, double x, double y )
{
    // This maps from a a logically rectangular 3D mesh to a sphere mesh using the mapping by:
    //    Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
    //       and Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains",
    //       SIAM REVIEW, Vol. 50, No. 4, pp.723–752 (2008)
    const double xc = 2 * x - 1; // Change domain to [-1,1]
    const double yc = 2 * y - 1; // Change domain to [-1,1]
    double xp, yp;
    if ( fabs( xc ) >= fabs( yc ) ) {
        auto p = map_c2p( method, fabs( xc ), fabs( yc ) );
        xp     = p.first;
        yp     = p.second;
    } else {
        auto p = map_c2p( method, fabs( yc ), fabs( xc ) );
        xp     = p.second;
        yp     = p.first;
    }
    if ( xc < 0.0 )
        xp = -xp;
    if ( yc < 0.0 )
        yp = -yp;
    return std::make_pair( r * xp, r * yp );
}
std::pair<double, double> map_circle_logical( double r, int method, double x, double y )
{
    // Get the points in the unit circle
    double xp = x / r;
    double yp = y / r;
    // Perform the inverse mapping to [-1,1]
    double xc = 0;
    double yc = 0;
    if ( fabs( xp ) >= fabs( yp ) ) {
        auto p = map_p2c( method, fabs( xp ), fabs( yp ) );
        xc     = p.first;
        yc     = p.second;
    } else {
        auto p = map_p2c( method, fabs( yp ), fabs( xp ) );
        xc     = p.second;
        yc     = p.first;
    }
    if ( xp < 0.0 )
        xc = -xc;
    if ( yp < 0.0 )
        yc = -yc;
    // Change domain to [0,1]
    return std::make_pair( 0.5 * ( xc + 1 ), 0.5 * ( yc + 1 ) );
}


/****************************************************************
 * Helper function to map x,y,z logical coordinates in [0,1]     *
 * to x,y,z coordinates in a sphere with radius r                *
 ****************************************************************/
AMP::Geometry::Point map_logical_sphere( double r, double x, double y, double z )
{
    constexpr double sqrt3 = 1.732050807568877; // sqrt(3)
    // This maps from a a logically rectangular 3D mesh to a sphere mesh using the mapping by:
    //    Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
    //       and Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains",
    //       SIAM REVIEW, Vol. 50, No. 4, pp.723–752 (2008)
    double xc = 2 * x - 1; // Change domain to [-1,1]
    double yc = 2 * y - 1; // Change domain to [-1,1]
    double zc = 2 * z - 1; // Change domain to [-1,1]
    double d  = std::max( { fabs( xc ), fabs( yc ), fabs( zc ) } );
    double r2 = sqrt( xc * xc + yc * yc + zc * zc );
    r2        = std::max( r2, 1e-10 );
    double d2 = d * d;
    double c  = r * ( d * d2 / r2 + ( 1.0 - d2 ) / sqrt3 );
    double x2 = c * xc;
    double y2 = c * yc;
    double z2 = c * zc;
    return { x2, y2, z2 };
}
AMP::Geometry::Point map_sphere_logical( double r, double x2, double y2, double z2 )
{
    constexpr double sqrt1_3 = 0.577350269189626; // 1/sqrt(3)
    // This maps from a physical coordinates to logical coordinates using the mapping by:
    //    Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
    //       and Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains",
    //       SIAM REVIEW, Vol. 50, No. 4, pp.723–752 (2008)
    double d1  = std::max( { fabs( x2 ), fabs( y2 ), fabs( z2 ) } );
    double r21 = sqrt( x2 * x2 + y2 * y2 + z2 * z2 );
    r21        = std::max( r21, 1e-10 );
    // Solve c = a + b/c^2
    double a  = r * sqrt1_3;
    double a2 = a * a;
    double b  = r * d1 * d1 * ( d1 / r21 - sqrt1_3 );
    double c  = 0;
    if ( fabs( b ) < 1e-8 * a ) {
        c = a;
        c = a + b / ( c * c );
    } else {
        c = 1.5 * sqrt( 12 * a2 * a * b + 81 * b * b ) + a2 * a + 13.5 * b;
        c = cbrt( c );
        c = 0.33333333333333333 * ( c + a2 / c + a );
    }
    // Compute (x,y,z)
    double xc = x2 / c;
    double yc = y2 / c;
    double zc = z2 / c;
    double x  = 0.5 + 0.5 * xc;
    double y  = 0.5 + 0.5 * yc;
    double z  = 0.5 + 0.5 * zc;
    return { x, y, z };
}


/****************************************************************
 * Helper function to map x,y logical coordinates in [0,1]       *
 * to x,y,z coordinates on the surface of a sphere               *
 ****************************************************************/
AMP::Geometry::Point map_logical_sphere_surface( double R, double x, double y )
{
    // This maps from a a logically rectangular 3D mesh to the surface of a sphere using:
    // Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids and
    //    Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains",
    //    SIAM REVIEW, Vol. 50, No. 4, pp. 723–752 (2008)
    double x2 = 2 * x - 1;  // Change domain to [-1,1]
    double x3 = fabs( x2 ); // We need to make x go from 1:0:1
    // Map x,y to the unit circle
    auto point = map_logical_circle( 1.0, 3, x3, y );
    double xp  = point.first;
    double yp  = point.second;
    double zp  = sqrt( fabs( 1.0 - ( xp * xp + yp * yp ) ) );
    if ( x2 < 0 )
        zp = -zp; // negate z in lower hemisphere
    xp *= R;
    yp *= R;
    zp *= R;
    return AMP::Geometry::Point( xp, yp, zp );
}
std::pair<double, double> map_sphere_surface_logical( double R, double x, double y, double z )
{
    double xp  = x / R;
    double yp  = y / R;
    auto point = map_circle_logical( 1.0, 3, xp, yp );
    double x2  = z < 0 ? -point.first : point.first;
    return { 0.5 + 0.5 * x2, point.second };
}


/****************************************************************
 * Helper function to map x,y,z logical coordinates in [0,1]     *
 * to x,y,z coordinates in a shell with r1 <= r <= r2            *
 ****************************************************************/
AMP::Geometry::Point map_logical_shell( double r1, double r2, double x, double y, double z )
{
    // This maps from a a logically rectangular 3D mesh to a shell mesh using the mapping by:
    // Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids and
    //    Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains",
    //    SIAM REVIEW, Vol. 50, No. 4, pp. 723–752 (2008)
    double Rz = r1 + z * ( r2 - r1 ); // radius based on z[0,1]
    return map_logical_sphere_surface( Rz, x, y );
}
AMP::Geometry::Point map_shell_logical( double r1, double r2, double x0, double y0, double z0 )
{
    // This maps from a a logically rectangular 3D mesh to a shell mesh using the mapping by:
    // Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids and
    //    Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains",
    //    SIAM REVIEW, Vol. 50, No. 4, pp. 723–752 (2008)
    double R = sqrt( x0 * x0 + y0 * y0 + z0 * z0 );
    auto xy  = map_sphere_surface_logical( R, x0, y0, z0 );
    double z = ( R - r1 ) / ( r2 - r1 );
    return { std::get<0>( xy ), std::get<1>( xy ), z };
}


/****************************************************************
 * Compute the distance to the surface of a plane                *
 ****************************************************************/
double distanceToPlane( const AMP::Mesh::Point &n,
                        const AMP::Mesh::Point &p0,
                        const AMP::Mesh::Point &pos,
                        const AMP::Mesh::Point &ang )
{
    double d = dot( n, ang );
    if ( fabs( d ) > 1e-6 ) {
        auto v = p0 - pos;
        auto t = dot( v, n ) / d;
        if ( t >= 0 )
            return t;
    }
    return std::numeric_limits<double>::infinity();
}


/****************************************************************
 * Compute the distance to a circle                              *
 ****************************************************************/
double distanceToCircle( double r, const AMP::Mesh::Point &pos, const AMP::Mesh::Point &ang )
{
    double r2   = pos.x() * pos.x() + pos.y() * pos.y();
    bool inside = r2 < r * r;
    double dx   = ang.x();
    double dy   = ang.y();
    double dr   = sqrt( dx * dx + dy * dy );
    double D    = pos.x() * ( pos.y() + ang.y() ) - ( pos.x() + ang.x() ) * pos.y();
    double t    = r * r * dr * dr - D * D;
    if ( t < 0 )
        return std::numeric_limits<double>::infinity();
    t = sqrt( t );
    double x1, x2, d1, d2;
    if ( ang.x() != 0.0 ) {
        double s = dy < 0 ? -1 : 1;
        x1       = ( D * dy + s * dx * t ) / ( dr * dr );
        x2       = ( D * dy - s * dx * t ) / ( dr * dr );
        d1       = ( x1 - pos.x() ) / ang.x();
        d2       = ( x2 - pos.x() ) / ang.x();
    } else {
        double s = dx < 0 ? -1 : 1;
        x1       = ( D * dx + s * dy * t ) / ( dr * dr );
        x2       = ( D * dx - s * dy * t ) / ( dr * dr );
        d1       = ( x1 - pos.y() ) / ang.y();
        d2       = ( x2 - pos.y() ) / ang.y();
    }
    if ( d1 < 0 )
        d1 = std::numeric_limits<double>::infinity();
    if ( d2 < 0 )
        d2 = std::numeric_limits<double>::infinity();
    double d = std::min( d1, d2 );
    return ( inside ? -1 : 1 ) * d;
}


/****************************************************************
 * Compute the distance to the surface of a cylinder             *
 ****************************************************************/
double distanceToCylinder( double r, double h, const Point &pos, const Point &ang )
{
    // Check if the point is inside the cylinder
    double r2   = pos.x() * pos.x() + pos.y() * pos.y();
    bool inside = std::abs( pos.z() ) <= 0.5 * h && r2 <= r * r;
    // First check the distance to the faces
    double d1 = ( 0.5 * h - pos.z() ) / ang.z();
    double d2 = ( -0.5 * h - pos.z() ) / ang.z();
    if ( d1 <= 0 )
        d1 = std::numeric_limits<double>::infinity();
    if ( d2 <= 0 )
        d2 = std::numeric_limits<double>::infinity();
    double d  = std::min( d1, d2 );
    auto pos2 = pos + d * ang;
    r2        = pos2.x() * pos2.x() + pos2.y() * pos2.y();
    if ( r2 <= r * r )
        return -d;
    else if ( ang.x() == 0 && ang.y() == 0 )
        return std::numeric_limits<double>::infinity();
    // Compute the intersection of a line with the circle of the cylinder
    d = std::abs( distanceToCircle( r, pos, ang ) );
    // Check that the z-point is within the cylinder
    double z = pos.z() + d * ang.z();
    if ( z < -0.5 * h || z > 0.5 * h )
        d = std::numeric_limits<double>::infinity();
    return ( inside ? -1 : 1 ) * d;
}


/****************************************************************
 * Compute the distance to the surface of a sphere               *
 *    http://paulbourke.net/geometry/circlesphere                *
 ****************************************************************/
double distanceToSphere( double r, const Point &pos, const Point &ang )
{
    double a = ang.x() * ang.x() + ang.y() * ang.y() + ang.z() * ang.z();
    double b = 2 * ( ang.x() * pos.x() + ang.y() * pos.y() + ang.z() * pos.z() );
    double c = pos.x() * pos.x() + pos.y() * pos.y() + pos.z() * pos.z() - r * r;
    double t = b * b - 4 * a * c;
    if ( t < 0 )
        return std::numeric_limits<double>::infinity();
    t         = sqrt( t );
    double d1 = ( -b + t ) / ( 2 * a );
    double d2 = ( -b - t ) / ( 2 * a );
    if ( d1 < 0 )
        d1 = std::numeric_limits<double>::infinity();
    if ( d2 < 0 )
        d2 = std::numeric_limits<double>::infinity();
    double d    = std::min( d1, d2 );
    double r2   = pos.x() * pos.x() + pos.y() * pos.y() + pos.z() * pos.z();
    bool inside = r2 < r * r;
    return ( inside ? -1 : 1 ) * d;
}


/****************************************************************
 * Compute the distance to the surface of a cone                *
 * http://lousodrome.net/blog/light/2017/01/03/intersection-of-a-ray-and-a-cone
 ****************************************************************/
double distanceToCone( const AMP::Mesh::Point &V,
                       double theta,
                       const AMP::Mesh::Point &pos,
                       const AMP::Mesh::Point &ang )
{
    double DV  = dot( ang, V );
    double CV  = dot( pos, V );
    double ct  = cos( theta );
    double ct2 = ct * ct;
    double a   = DV * DV - ct2;
    double b   = 2 * ( DV * CV - dot( ang, pos ) * ct2 );
    double c   = CV * CV - dot( pos, pos ) * ct2;
    double t   = b * b - 4 * a * c;
    if ( t < 0 )
        return std::numeric_limits<double>::infinity();
    t         = sqrt( t );
    double d1 = ( -b + t ) / ( 2 * a );
    double d2 = ( -b - t ) / ( 2 * a );
    if ( d1 < 0 )
        d1 = std::numeric_limits<double>::infinity();
    if ( d2 < 0 )
        d2 = std::numeric_limits<double>::infinity();
    if ( dot( pos + d1 * ang, V ) < 0 )
        d1 = std::numeric_limits<double>::infinity();
    if ( dot( pos + d2 * ang, V ) < 0 )
        d2 = std::numeric_limits<double>::infinity();
    return std::min( d1, d2 );
}


} // namespace GeometryHelpers
} // namespace Geometry
} // namespace AMP
