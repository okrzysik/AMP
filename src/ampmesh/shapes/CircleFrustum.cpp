#include "AMP/ampmesh/shapes/CircleFrustum.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>
#include <vector>


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructors                                          *
 ********************************************************/
CircleFrustum::CircleFrustum( std::shared_ptr<AMP::Database> db ) : d_offset{ 0, 0, 0 }
{
    d_r[0]   = db->getScalar<double>( "BaseRadius" );
    d_r[1]   = db->getScalar<double>( "TopRadius" );
    d_h      = db->getScalar<double>( "Height" );
    auto dir = db->getString( "Dir" );
    if ( dir == "-x" )
        d_dir = 0;
    else if ( dir == "+x" )
        d_dir = 1;
    else if ( dir == "-y" )
        d_dir = 2;
    else if ( dir == "+y" )
        d_dir = 3;
    else if ( dir == "-z" )
        d_dir = 4;
    else if ( dir == "+z" )
        d_dir = 5;
    else
        AMP_ERROR( "Invalid value for Dir" );
    d_physicalDim = 3;
    d_logicalDim  = 3;
    AMP_INSIST( d_r[0] > d_r[1] && d_r[1] > 0, "Invalid value for r" );
    // Compute the apex of the underlying cone
    double h2      = d_h * d_r[0] / ( d_r[0] - d_r[1] );
    d_C            = { 0, 0, 0 };
    d_C[d_dir / 2] = d_dir % 2 == 0 ? -h2 : h2;
    d_theta        = atan( d_r[0] / h2 );
}
CircleFrustum::CircleFrustum( const std::array<double, 2> &r, int dir, double height )
    : LogicalGeometry(), d_dir( dir ), d_r{ r[0], r[1] }, d_h( height ), d_offset{ 0, 0, 0 }
{
    d_physicalDim = 3;
    d_logicalDim  = 3;
    AMP_INSIST( d_r[0] > d_r[1] && d_r[1] > 0, "Invalid value for r" );
    AMP_INSIST( d_dir < 6, "Invalid value for dir" );
    // Compute the apex of the underlying cone
    double h2      = d_h * d_r[0] / ( d_r[0] - d_r[1] );
    d_C            = { 0, 0, 0 };
    d_C[d_dir / 2] = d_dir % 2 == 0 ? -h2 : h2;
    d_theta        = atan( d_r[0] / h2 );
}


/********************************************************
 * Compute the nearest point on the surface              *
 ********************************************************/
Point CircleFrustum::nearest( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return {};
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
    double d1 = -p0[dir2] / ang[dir2];
    double d2 = ( s * d_h - p0[dir2] ) / ang[dir2];
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
        // We are at the - face
        return 0;
    } else if ( std::abs( 1 - p.z() ) < t ) {
        // We are at the + face
        return 1;
    } else {
        // We are at the cone face
        return 2;
    }
}
Point CircleFrustum::surfaceNorm( const Point &pos ) const
{
    int s    = surface( pos );
    Point v  = { 0, 0, 0 };
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2];
    if ( s == 0 ) {
        v[d_dir / 2] = d_dir % 2 == 0 ? 1 : -1;
    } else if ( s == 1 ) {
        v[d_dir / 2] = d_dir % 2 == 0 ? -1 : 1;
    } else {
        double sin_t = sin( d_theta );
        double cos_t = cos( d_theta );
        if ( d_dir == 0 ) {
            double r = sqrt( y * y + z * z );
            v        = { -sin_t, cos_t * y / r, cos_t * z / r };
        } else if ( d_dir == 1 ) {
            double r = sqrt( y * y + z * z );
            v        = { sin_t, cos_t * y / r, cos_t * z / r };
        } else if ( d_dir == 2 ) {
            double r = sqrt( x * x + z * z );
            v        = { cos_t * x / r, -sin_t, cos_t * z / r };
        } else if ( d_dir == 3 ) {
            double r = sqrt( x * x + z * z );
            v        = { cos_t * x / r, sin_t, cos_t * z / r };
        } else if ( d_dir == 4 ) {
            double r = sqrt( x * x + y * y );
            v        = { cos_t * x / r, cos_t * y / r, -sin_t };
        } else {
            double r = sqrt( x * x + y * y );
            v        = { cos_t * x / r, cos_t * y / r, sin_t };
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
 * Return the logical grid                               *
 ********************************************************/
std::vector<int> CircleFrustum::getLogicalGridSize( const std::vector<int> &x ) const
{
    AMP_INSIST( x.size() == 2, "Size must be an array of length 2" );
    return { x[0], x[0], x[1] };
}
std::vector<int> CircleFrustum::getLogicalGridSize( const std::vector<double> &res ) const
{
    AMP_INSIST( res.size() == 3u, "Resolution must be an array of length 3" );
    AMP_ERROR( "Not finished" );
    return {};
}
std::vector<bool> CircleFrustum::getPeriodicDim() const { return { false, false, false }; }
std::vector<int> CircleFrustum::getLogicalSurfaceIds() const { return { 2, 2, 2, 2, 0, 1 }; }


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void CircleFrustum::displace( const double *x )
{
    d_offset[0] += x[0];
    d_offset[1] += x[1];
    d_offset[2] += x[2];
}


/********************************************************
 * Clone the object                                      *
 ********************************************************/
std::unique_ptr<AMP::Geometry::Geometry> CircleFrustum::clone() const
{
    return std::make_unique<CircleFrustum>( *this );
}

} // namespace Geometry
} // namespace AMP
