#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/ampmesh/Geometry.h"
#include "AMP/utils/Utilities.h"
#if defined( USING_ICC )
#include "AMP/utils/UtilityMacros.h"
#endif

#include <algorithm>
#include <cmath>


#if defined( USING_ICC )
DISABLE_WARNINGS
#endif

/****************************************************************
 * Overload basic array operations                               *
 ****************************************************************/
template<size_t N>
static inline std::array<double, N> operator*( double x, const std::array<double, N> &y )
{
    if constexpr ( N == 1 )
        return { x * y[0] };
    else if constexpr ( N == 2 )
        return { x * y[0], x * y[1] };
    else if constexpr ( N == 3 )
        return { x * y[0], x * y[1], x * y[2] };
}
template<size_t N>
static inline std::array<double, N> operator+( const std::array<double, N> &x,
                                               const std::array<double, N> &y )
{
    if constexpr ( N == 1 )
        return { x[0] + y[0] };
    else if constexpr ( N == 2 )
        return { x[0] + y[0], x[1] + y[1] };
    else if constexpr ( N == 3 )
        return { x[0] + y[0], x[1] + y[1], x[2] + y[2] };
}
template<size_t N>
static inline std::array<double, N> operator-( const std::array<double, N> &x,
                                               const std::array<double, N> &y )
{
    if constexpr ( N == 1 )
        return { x[0] - y[0] };
    else if constexpr ( N == 2 )
        return { x[0] - y[0], x[1] - y[1] };
    else if constexpr ( N == 3 )
        return { x[0] - y[0], x[1] - y[1], x[2] - y[2] };
}
template<size_t N>
static inline double abs( const std::array<double, N> &x )
{
    if constexpr ( N == 1 )
        return std::abs( x[0] );
    else if constexpr ( N == 2 )
        return sqrt( x[0] * x[0] + x[1] * x[1] );
    else if constexpr ( N == 3 )
        return sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
}
template<size_t N>
static inline double dot( const std::array<double, N> &x, const std::array<double, N> &y )
{
    if constexpr ( N == 1 )
        return x[0] * y[0];
    else if constexpr ( N == 2 )
        return x[0] * y[0] + x[1] * y[1];
    else if constexpr ( N == 3 )
        return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}
static inline std::array<double, 3> cross( const std::array<double, 3> &x,
                                           const std::array<double, 3> &y )
{
    return { x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0] };
}
template<size_t N>
static inline std::array<double, N> normalize( const std::array<double, N> &x )
{
    if constexpr ( N == 1 ) {
        return { 1.0 };
    } else if constexpr ( N == 2 ) {
        double tmp = 1.0 / sqrt( x[0] * x[0] + x[1] * x[1] );
        return { tmp * x[0], tmp * x[1] };
    } else if constexpr ( N == 3 ) {
        double tmp = 1.0 / sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
        return { tmp * x[0], tmp * x[1], tmp * x[2] };
    }
}


namespace AMP::Geometry::GeometryHelpers {


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
 * Compute the distance to a box                                 *
 ****************************************************************/
double distanceToBox( const AMP::Mesh::Point &pos,
                      const AMP::Mesh::Point &ang,
                      const std::array<double, 6> &range )
{
    constexpr double tol = 1e-12;
    // Compute the distance to each surface
    double d1 = ( range[0] - pos.x() ) / ang.x();
    double d2 = ( range[1] - pos.x() ) / ang.x();
    double d3 = ( range[2] - pos.y() ) / ang.y();
    double d4 = ( range[3] - pos.y() ) / ang.y();
    double d5 = ( range[4] - pos.z() ) / ang.z();
    double d6 = ( range[5] - pos.z() ) / ang.z();
    if ( d1 < 0 )
        d1 = std::numeric_limits<double>::infinity();
    if ( d2 < 0 )
        d2 = std::numeric_limits<double>::infinity();
    if ( d3 < 0 )
        d3 = std::numeric_limits<double>::infinity();
    if ( d4 < 0 )
        d4 = std::numeric_limits<double>::infinity();
    if ( d5 < 0 )
        d5 = std::numeric_limits<double>::infinity();
    if ( d6 < 0 )
        d6 = std::numeric_limits<double>::infinity();
    // Check if the intersection of each surface is within the bounds of the box
    auto p1     = pos + d1 * ang;
    auto p2     = pos + d2 * ang;
    auto p3     = pos + d3 * ang;
    auto p4     = pos + d4 * ang;
    auto p5     = pos + d5 * ang;
    auto p6     = pos + d6 * ang;
    auto inside = [range]( const Point &p ) {
        return ( ( p.x() >= range[0] - tol ) && ( p.x() <= range[1] + tol ) ) &&
               ( ( p.y() >= range[2] - tol ) && ( p.y() <= range[3] + tol ) ) &&
               ( ( p.z() >= range[4] - tol ) && ( p.z() <= range[5] + tol ) );
    };
    if ( !inside( p1 ) )
        d1 = std::numeric_limits<double>::infinity();
    if ( !inside( p2 ) )
        d2 = std::numeric_limits<double>::infinity();
    if ( !inside( p3 ) )
        d3 = std::numeric_limits<double>::infinity();
    if ( !inside( p4 ) )
        d4 = std::numeric_limits<double>::infinity();
    if ( !inside( p5 ) )
        d5 = std::numeric_limits<double>::infinity();
    if ( !inside( p6 ) )
        d6 = std::numeric_limits<double>::infinity();
    // Return the closest surface
    double d = std::min( { d1, d2, d3, d4, d5, d6 } );
    if ( inside( pos ) && d < 1e100 )
        d = -d;
    return d;
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
    if ( d1 <= 0 || ang.z() == 0 )
        d1 = std::numeric_limits<double>::infinity();
    if ( d2 <= 0 || ang.z() == 0 )
        d2 = std::numeric_limits<double>::infinity();
    auto pos2 = pos + std::min( d1, d2 ) * ang;
    r2        = pos2.x() * pos2.x() + pos2.y() * pos2.y();
    if ( r2 > r * r ) {
        d1 = std::numeric_limits<double>::infinity();
        d2 = std::numeric_limits<double>::infinity();
    }
    // Compute the intersection of a line with the circle of the cylinder
    double d3 = std::abs( distanceToCircle( r, pos, ang ) );
    // Check that the z-point is within the cylinder
    double z = pos.z() + d3 * ang.z();
    if ( fabs( z ) > 0.5 * h )
        d3 = std::numeric_limits<double>::infinity();
    // Return the closest surface
    double d = std::min( { d1, d2, d3 } );
    if ( d < 1e100 )
        return ( inside ? -1 : 1 ) * d;
    return std::numeric_limits<double>::infinity();
}


/****************************************************************
 * Compute the distance to the surface of a cylinder             *
 ****************************************************************/
double distanceToTube( double r_min, double r_max, double h, const Point &pos, const Point &ang )
{
    double r_min2 = r_min * r_min;
    double r_max2 = r_max * r_max;
    // Check if the point is inside the tube
    double r2   = pos.x() * pos.x() + pos.y() * pos.y();
    bool inside = std::abs( pos.z() ) <= 0.5 * h && r2 >= r_min2 && r2 <= r_max2;
    // Check the distance to the faces
    double d1 = ( 0.5 * h - pos.z() ) / ang.z();
    double d2 = ( -0.5 * h - pos.z() ) / ang.z();
    if ( d1 <= 0 || ang.z() == 0 )
        d1 = std::numeric_limits<double>::infinity();
    if ( d2 <= 0 || ang.z() == 0 )
        d2 = std::numeric_limits<double>::infinity();
    auto pos2 = pos + std::min( d1, d2 ) * ang;
    r2        = pos2.x() * pos2.x() + pos2.y() * pos2.y();
    if ( r2 < r_min2 || r2 > r_max2 ) {
        d1 = std::numeric_limits<double>::infinity();
        d2 = std::numeric_limits<double>::infinity();
    }
    // Check the intersection of a line with the circles of the tube
    auto checkCircle = [pos, ang, h]( double r ) {
        if ( ang.x() == 0 && ang.y() == 0 )
            return std::numeric_limits<double>::infinity();
        double d = std::abs( distanceToCircle( r, pos, ang ) );
        double z = pos.z() + d * ang.z();
        // We did not intersect with the surface, check for a second intersection
        if ( fabs( z ) > 0.5 * h ) {
            d += 1e-8;
            auto pos2 = pos + d * ang;
            double d2 = std::abs( distanceToCircle( r, pos2, ang ) );
            d += d2;
            z = pos.z() + d * ang.z();
            if ( fabs( z ) > 0.5 * h )
                d = std::numeric_limits<double>::infinity();
        }
        return d;
    };
    double d3 = checkCircle( r_min );
    double d4 = checkCircle( r_max );
    // Return the closest surface
    double d = std::min( { d1, d2, d3, d4 } );
    if ( d < 1e100 )
        return ( inside ? -1 : 1 ) * d;
    return std::numeric_limits<double>::infinity();
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


/****************************************************************
 * Compute the normal to a plane                                 *
 ****************************************************************/
std::array<double, 3> normal( const std::array<double, 3> &v1,
                              const std::array<double, 3> &v2,
                              const std::array<double, 3> &v3 )
{
    return normalize( cross( v1 - v2, v1 - v3 ) );
}


/****************************************************************
 * Get the Barycentric coordinates (u, v, w) for point p with   *
 *   respect to triangle (a, b, c)                              *
 ***************************************************************/
template<>
std::array<double, 3> barycentric<3, 3>( const std::array<double, 3> ( &x )[3],
                                         const std::array<double, 3> &p )
{
    auto v0  = x[1] - x[0];
    auto v1  = x[2] - x[0];
    auto v2  = p - x[0];
    auto d00 = dot( v0, v0 );
    auto d01 = dot( v0, v1 );
    auto d11 = dot( v1, v1 );
    auto d20 = dot( v2, v0 );
    auto d21 = dot( v2, v1 );
    auto inv = 1.0 / ( d00 * d11 - d01 * d01 );
    auto v   = inv * ( d11 * d20 - d01 * d21 );
    auto w   = inv * ( d00 * d21 - d01 * d20 );
    auto u   = 1.0 - v - w;
    return { u, v, w };
}


/****************************************************************
 * Compute the nearest point to a line segment                   *
 ****************************************************************/
std::array<double, 3> nearest( const std::array<double, 3> &A,
                               const std::array<double, 3> &B,
                               const std::array<double, 3> &P )
{
    auto v = B - A;
    auto u = A - P;
    auto t = -dot( v, u ) / dot( v, v );
    t      = std::min( std::max( t, 0.0 ), 1.0 );
    return ( 1.0 - t ) * A + t * B;
}


/****************************************************************
 * Compute the nearest point to a triangle                       *
 ****************************************************************/
std::array<double, 3> nearest( const std::array<double, 3> ( &v )[3],
                               const std::array<double, 3> &p0 )
{
    constexpr double TOL = 1e-8;
    // Get the normal and a point on the plane containing the triangle
    auto n = AMP::Geometry::GeometryHelpers::normal( v[0], v[1], v[2] );
    // Find the closest point on the plane
    double d = dot( n, p0 - v[0] );
    auto p   = p0 - d * n;
    // Compute the barycentric coordinates
    auto L = barycentric<3, 3>( v, p );
    if ( L[0] > -TOL && L[1] > -TOL && L[2] > -TOL ) {
        // Point is inside triangle
        return p;
    } else if ( L[1] <= 0 && L[2] <= 0 ) {
        // Point is closest to the first vertex
        return v[0];
    } else if ( L[0] <= 0 && L[2] <= 0 ) {
        // Point is closest to the second vertex
        return v[1];
    } else if ( L[0] <= 0 && L[1] <= 0 ) {
        // Point is closest to the third vertex
        return v[2];
    } else if ( L[0] <= 0 ) {
        // Point is closest to line between second and third verticies
        return nearest( v[1], v[2], p0 );
    } else if ( L[1] <= 0 ) {
        // Point is closest to line between first and third verticies
        return nearest( v[0], v[2], p0 );
    } else {
        // Point is closest to line between first and second verticies
        return nearest( v[0], v[1], p0 );
    }
}


/****************************************************************
 * Recursively divide the triangle                               *
 ****************************************************************/
std::vector<AMP::Mesh::Point> subdivide( const std::array<AMP::Mesh::Point, 3> &v, double res )
{
    double d1  = ( v[0] - v[1] ).norm();
    double d2  = ( v[0] - v[2] ).norm();
    double d3  = ( v[1] - v[2] ).norm();
    double max = std::max( { d1, d2, d3 } );
    if ( max < res * res )
        return std::vector<AMP::Mesh::Point>();
    std::vector<AMP::Mesh::Point> s1, s2;
    if ( d1 == max ) {
        auto p = 0.5 * ( v[0] + v[1] );
        s1     = subdivide( { v[0], p, v[2] }, res );
        s2     = subdivide( { p, v[1], v[2] }, res );
        s1.push_back( p );
    } else if ( d2 == max ) {
        auto p = 0.5 * ( v[0] + v[2] );
        s1     = subdivide( { v[0], v[1], p }, res );
        s2     = subdivide( { p, v[1], v[2] }, res );
        s1.push_back( p );
    } else {
        auto p = 0.5 * ( v[1] + v[2] );
        s1     = subdivide( { v[0], v[1], p }, res );
        s2     = subdivide( { v[0], p, v[2] }, res );
        s1.push_back( p );
    }
    const double tol = 0.5 * res * res;
    for ( const auto &p : s2 ) {
        bool found = false;
        for ( const auto &p0 : s1 )
            found = found || ( p0 - p ).norm() < tol;
        if ( !found )
            s1.push_back( p );
    }
    return s1;
}


/****************************************************************
 * Compute the normal to the plane formed by 3 points            *
 ****************************************************************/
AMP::Mesh::Point
normal( const AMP::Mesh::Point &a, const AMP::Mesh::Point &b, const AMP::Mesh::Point &c )
{
    return normalize( cross( b - a, c - a ) );
}


} // namespace AMP::Geometry::GeometryHelpers

#if defined( USING_ICC )
ENABLE_WARNINGS
#endif
