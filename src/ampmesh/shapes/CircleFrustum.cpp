#include "AMP/ampmesh/shapes/CircleFrustum.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>
#include <vector>


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructors                                          *
 ********************************************************/
CircleFrustum::CircleFrustum( const std::array<double, 2> &r, int dir, double height )
    : Geometry(), d_dir( dir ), d_r{ r[0], r[1] }, d_h( height ), d_offset{ 0, 0, 0 }
{
    d_physicalDim = 3;
    d_logicalDim  = 3;
    AMP_INSIST( r[0] > r[1] && r[1] > 0, "Invalid value for r" );
    AMP_INSIST( dir >= 0 && dir < 6, "Invalid value for dir" );
    // Compute the apex of the underlying cone
    double h2    = height * r[0] / ( r[0] - r[1] );
    d_C          = { 0, 0, 0 };
    d_C[dir / 2] = d_dir % 2 == 0 ? -h2 : h2;
    d_theta      = atan( r[0] / h2 );
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double CircleFrustum::distance( const Point &pos, const Point &ang ) const
{
    auto dir2 = d_dir / 2;
    // Remove the offset
    Point p0 = pos;
    p0.x() -= d_offset[0];
    p0.y() -= d_offset[1];
    p0.z() -= d_offset[2];
    // Compute the intersection with the infinite cone
    Point V  = { 0, 0, 0 };
    V[dir2]  = d_dir % 2 == 0 ? 1 : -1;
    double d = std::abs( GeometryHelpers::distanceToCone( V, d_theta, p0 - d_C, ang ) );
    auto p   = p0 + d * ang;
    if ( d_dir == 0 && ( p.x() < -d_h || p.x() > 0 ) ) {
        d = std::numeric_limits<double>::infinity();
    } else if ( d_dir == 1 && ( p.x() < 0 || p.x() > d_h ) ) {
        d = std::numeric_limits<double>::infinity();
    } else if ( d_dir == 2 && ( p.y() < -d_h || p.y() > 0 ) ) {
        d = std::numeric_limits<double>::infinity();
    } else if ( d_dir == 3 && ( p.y() < 0 || p.y() > d_h ) ) {
        d = std::numeric_limits<double>::infinity();
    } else if ( d_dir == 4 && ( p.z() < -d_h || p.z() > 0 ) ) {
        d = std::numeric_limits<double>::infinity();
    } else if ( d_dir == 5 && ( p.z() < 0 || p.z() > d_h ) ) {
        d = std::numeric_limits<double>::infinity();
    }
    // Compute the intersection with the planes slicing the cone
    bool swap = d_dir == 0 || d_dir == 1 || d_dir == 3 || d_dir == 4;
    double s  = swap ? -1 : 1;
    double d1 = -p0[dir2] / ( s * ang[dir2] );
    double d2 = ( s * d_h - p0[dir2] ) / ( s * ang[dir2] );
    auto p1   = p0 + d1 * ang;
    auto p2   = p0 + d2 * ang;
    double r1, r2;
    if ( dir2 == 0 ) {
        r1 = p1.y() * p1.y() + p1.z() * p1.z();
        r2 = p2.y() * p2.y() + p2.z() * p2.z();
    } else if ( dir2 == 1 ) {
        r1 = p1.x() * p1.x() + p1.z() * p1.z();
        r2 = p2.x() * p2.x() + p2.z() * p2.z();
    } else {
        r1 = p1.x() * p1.x() + p1.y() * p1.y();
        r2 = p2.x() * p2.x() + p2.y() * p2.y();
    }
    if ( ( d1 < 0 ) || ( r1 > d_r[0] * d_r[0] ) )
        d1 = std::numeric_limits<double>::infinity();
    if ( ( d2 < 0 ) || ( r2 > d_r[1] * d_r[1] ) )
        d2 = std::numeric_limits<double>::infinity();
    // Keep the closest intersection
    d = std::min( { d, d1, d2 } );
    // Check if the point is inside the volume
    if ( d < 1e200 ) {
        if ( inside( pos ) )
            d = -d;
    }
    return d;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool CircleFrustum::inside( const Point &pos ) const
{
    constexpr double t1 = -1e-12;
    constexpr double t2 = 1.0 + 1e-12;
    // Get and check the logical coordinates
    auto p = logical( pos );
    return p.x() >= t1 && p.x() <= t2 && p.y() >= t1 && p.y() <= t2 && p.z() >= t1 && p.z() <= t2;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int CircleFrustum::surface( const Point &pos ) const
{
    auto p   = logical( pos );
    double t = std::min(
        { std::abs( p.x() ), std::abs( 1 - p.x() ), std::abs( p.y() ), std::abs( 1 - p.y() ) } );
    if ( std::abs( p.z() ) < t ) {
        // We are at the large face
        bool swap = d_dir == 0 || d_dir == 1 || d_dir == 3 || d_dir == 4;
        return swap ? 1 : 0;
    } else if ( std::abs( 1 - p.z() ) < t ) {
        // We are at the small face
        bool swap = d_dir == 0 || d_dir == 1 || d_dir == 3 || d_dir == 4;
        return swap ? 0 : 1;
    } else {
        // We are at the cone face
        return 2;
    }
}
Point CircleFrustum::surfaceNorm( const Point &pos ) const
{
    int s   = surface( pos );
    Point v = { 0, 0, 0 };
    if ( s == 0 ) {
        v[d_dir / 2] = d_dir % 2 == 0 ? 1 : -1;
    } else if ( s == 1 ) {
        v[d_dir / 2] = d_dir % 2 == 0 ? -1 : 1;
    } else {
        if ( d_dir == 0 ) {
            double r = pos.y() * pos.y() + pos.z() * pos.z();
            v = { sin( d_theta ), cos( d_theta ) * pos.y() / r, cos( d_theta ) * pos.z() / r };
        } else if ( d_dir == 1 ) {
            double r = pos.y() * pos.y() + pos.z() * pos.z();
            v = { sin( d_theta ), cos( d_theta ) * pos.y() / r, cos( d_theta ) * pos.z() / r };
        } else if ( d_dir == 2 ) {
            double r = pos.x() * pos.x() + pos.z() * pos.z();
            v = { cos( d_theta ) * pos.x() / r, sin( d_theta ), cos( d_theta ) * pos.z() / r };
        } else if ( d_dir == 3 ) {
            double r = pos.x() * pos.x() + pos.z() * pos.z();
            v = { cos( d_theta ) * pos.x() / r, -sin( d_theta ), cos( d_theta ) * pos.z() / r };
        } else if ( d_dir == 4 ) {
            double r = pos.x() * pos.x() + pos.y() * pos.y();
            v = { cos( d_theta ) * pos.x() / r, cos( d_theta ) * pos.y() / r, sin( d_theta ) };
        } else {
            double r = pos.x() * pos.x() + pos.y() * pos.y();
            v = { cos( d_theta ) * pos.x() / r, cos( d_theta ) * pos.y() / r, -sin( d_theta ) };
        }
    }
    return v;
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point CircleFrustum::physical( const Point &pos ) const
{
    // Swap directions to preserve positive volume
    auto p0   = pos;
    bool swap = d_dir == 0 || d_dir == 1 || d_dir == 3 || d_dir == 4;
    if ( swap )
        p0.z() = 1.0 - p0.z();
    // Get the height and current radius
    Point p  = { 0, 0, 0 };
    p.z()    = d_h * p0.z();
    double r = d_r[0] * ( 1.0 - p0.z() ) + d_r[1] * p0.z();
    // Map the x/y coordinates to a circle
    std::tie( p.x(), p.y() ) = GeometryHelpers::map_logical_circle( r, 2, p0.x(), p0.y() );
    // Rotate the coordinates
    if ( d_dir == 0 ) {
        p = { -p.z(), p.x(), p.y() };
    } else if ( d_dir == 1 ) {
        p = { -p.z(), p.x(), p.y() };
    } else if ( d_dir == 2 ) {
        p = { p.x(), -p.z(), p.y() };
    } else if ( d_dir == 3 ) {
        p = { p.x(), p.z(), p.y() };
    } else if ( d_dir == 4 ) {
        p = { p.x(), p.y(), -p.z() };
    } else {
        p = { p.x(), p.y(), p.z() };
    }
    // Add the offset
    p.x() += d_offset[0];
    p.y() += d_offset[1];
    p.z() += d_offset[2];
    return p;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point CircleFrustum::logical( const Point &pos ) const
{
    // Remove the offset
    Point p0 = pos;
    p0.x() -= d_offset[0];
    p0.y() -= d_offset[1];
    p0.z() -= d_offset[2];
    // Get the physical point in a non-rotated frame
    if ( d_dir == 0 ) {
        p0 = { p0.y(), p0.z(), -p0.x() };
    } else if ( d_dir == 1 ) {
        p0 = { p0.y(), p0.z(), p0.x() };
    } else if ( d_dir == 2 ) {
        p0 = { p0.x(), p0.z(), -p0.y() };
    } else if ( d_dir == 3 ) {
        p0 = { p0.x(), p0.z(), p0.y() };
    } else if ( d_dir == 4 ) {
        p0 = { p0.x(), p0.y(), -p0.z() };
    } else {
        p0 = { p0.x(), p0.y(), p0.z() };
    }
    // Get the logical height and current radius
    Point p  = { 0, 0, 0 };
    p.z()    = p0.z() / d_h;
    double r = d_r[0] * ( 1.0 - p.z() ) + d_r[1] * p.z();
    // Map the x/y coordinates from a circle
    std::tie( p.x(), p.y() ) = GeometryHelpers::map_circle_logical( r, 2, p0.x(), p0.y() );
    // Swap directions to preserve positive volume
    bool swap = d_dir == 0 || d_dir == 1 || d_dir == 3 || d_dir == 4;
    if ( swap )
        p.z() = 1.0 - p.z();
    return p;
}


/********************************************************
 * Return the centroid and bounding box                  *
 ********************************************************/
Point CircleFrustum::centroid() const
{
    auto dir2 = d_dir / 2;
    Point p   = { 0, 0, 0 };
    p[dir2]   = ( d_dir % 2 == 0 ? -1 : 1 ) * d_h / 2;
    p.x() += d_offset[0];
    p.y() += d_offset[1];
    p.z() += d_offset[2];
    return p;
}
std::pair<Point, Point> CircleFrustum::box() const
{
    auto dir2 = d_dir / 2;
    Point lb  = { -d_r[0], -d_r[0], -d_r[0] };
    Point ub  = { d_r[0], d_r[0], d_r[0] };
    lb[dir2]  = d_dir % 2 == 0 ? -d_h : 0;
    ub[dir2]  = d_dir % 2 == 0 ? 0 : d_h;
    lb.x() += d_offset[0];
    lb.y() += d_offset[1];
    lb.z() += d_offset[2];
    ub.x() += d_offset[0];
    ub.y() += d_offset[1];
    ub.z() += d_offset[2];
    return { lb, ub };
}


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void CircleFrustum::displaceMesh( const double *x )
{
    d_offset[0] += x[0];
    d_offset[1] += x[1];
    d_offset[2] += x[2];
}


} // namespace Geometry
} // namespace AMP
