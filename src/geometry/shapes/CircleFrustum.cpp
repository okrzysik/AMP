#include "AMP/geometry/shapes/CircleFrustum.h"
#include "AMP/IO/HDF5.h"
#include "AMP/geometry/GeometryHelpers.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UtilityMacros.h"

#include <algorithm>
#include <vector>


namespace AMP::Geometry {


/********************************************************
 * Constructors                                          *
 ********************************************************/
CircleFrustum::CircleFrustum( std::shared_ptr<const AMP::Database> db )
    : d_dir( 0 ), d_r{ 0, 0 }, d_h( 0 ), d_offset{ 0, 0, 0 }
{
    double r1 = db->getScalar<double>( "BaseRadius" );
    double r2 = db->getScalar<double>( "TopRadius" );
    double h  = db->getScalar<double>( "Height" );
    auto dir  = db->getString( "Dir" );
    int dir2  = 0;
    if ( dir == "-x" )
        dir2 = 0;
    else if ( dir == "+x" )
        dir2 = 1;
    else if ( dir == "-y" )
        dir2 = 2;
    else if ( dir == "+y" )
        dir2 = 3;
    else if ( dir == "-z" )
        dir2 = 4;
    else if ( dir == "+z" )
        dir2 = 5;
    else
        AMP_ERROR( "Invalid value for Dir" );
    initialize( dir2, { r1, r2 }, h );
}
CircleFrustum::CircleFrustum( const std::array<double, 2> &r, int dir, double height )
    : LogicalGeometry(), d_dir( 0 ), d_r{ 0, 0 }, d_h( 0 ), d_offset{ 0, 0, 0 }
{
    initialize( dir, r, height );
}
void CircleFrustum::initialize( int dir, const std::array<double, 2> &r, double h )
{
    d_ids         = { 2, 2, 2, 2, 0, 1 };
    d_isPeriodic  = { false, false, false };
    d_dir         = dir;
    d_r[0]        = r[0];
    d_r[1]        = r[1];
    d_h           = h;
    d_offset[0]   = 0;
    d_offset[1]   = 0;
    d_offset[2]   = 0;
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
 * Note: this is currently an approximation              *
 ********************************************************/
Point CircleFrustum::nearest( const Point &pos ) const
{
#if 1
    auto L = logical( pos );
    L.x()  = std::max( L.x(), 0.0 );
    L.x()  = std::min( L.x(), 1.0 );
    L.y()  = std::max( L.y(), 0.0 );
    L.y()  = std::min( L.y(), 1.0 );
    L.z()  = std::max( L.z(), 0.0 );
    L.z()  = std::min( L.z(), 1.0 );
#else
    // Get the logical coordinates
    auto L   = logical( pos );
    auto tmp = physical( L );
    if ( L.z() < 0.0 || L.z() > 1.0 ) {
        // Get the radius and z-position
        double z  = d_h * L.z();
        double rz = d_r[0] * ( 1.0 - z ) + d_r[1] * z;
        auto pl   = GeometryHelpers::map_logical_circle( rz, 2, L.x(), L.y() );
        double r  = sqrt( pl[0] * pl[0] + pl[1] * pl[1] );
        if ( L.z() < 0.0 && r < d_r[0] ) {
            // Closest point is on the bottom surface
            L.z() = 0.0;
        } else if ( L.z() > 1.0 && r < d_r[1] ) {
            // Closest point is on the top surface
            L.z() = 1.0;
        } else {
            // Closest point is on the outer surface
            // Determine the closest point on the line of the radius
            std::array<double, 2> A = { 0.5 * d_r[0], 0.0 };
            std::array<double, 2> B = { 0.5 * d_r[1], d_h };
            auto p                  = GeometryHelpers::nearest( A, B, { r, z } );
            L.z()                   = p[1] / d_h;
        }
    }
    L.x() = std::max( L.x(), 0.0 );
    L.x() = std::min( L.x(), 1.0 );
    L.y() = std::max( L.y(), 0.0 );
    L.y() = std::min( L.y(), 1.0 );
#endif
    return physical( L );
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
    double r = d_r[0] * ( 1.0 - p0.z() ) + d_r[1] * p0.z();
    double z = d_h * p0.z();
    // Map the x/y coordinates to a circle
    auto pl = GeometryHelpers::map_logical_circle( r, 2, p0.x(), p0.y() );
    Point p = { pl[0], pl[1], z };
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
    double z = p0.z() / d_h;
    double r = d_r[0] * ( 1.0 - z ) + d_r[1] * z;
    // Map the x/y coordinates from a circle
    auto pl = GeometryHelpers::map_circle_logical( r, 2, p0.x(), p0.y() );
    Point p = { pl[0], pl[1], z };
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
    double R11 = d_r[0] * d_r[0];
    double R12 = d_r[0] * d_r[1];
    double R22 = d_r[1] * d_r[1];
    double z   = 0.25 * d_h * ( R11 + 2 * R12 + 3 * R22 ) / ( R11 + R12 + R22 );
    auto dir2  = d_dir / 2;
    Point p    = { 0, 0, 0 };
    p[dir2]    = ( d_dir % 2 == 0 ? -1 : 1 ) * z;
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
 * Return the volume                                     *
 ********************************************************/
double CircleFrustum::volume() const
{
    constexpr double pi = 3.141592653589793;
    double R            = 0.5 * ( d_r[0] + d_r[1] );
    return d_h * pi * R * R;
}


/********************************************************
 * Return the logical grid                               *
 ********************************************************/
std::vector<int> CircleFrustum::getLogicalGridSize( const std::vector<int> &x ) const
{
    AMP_INSIST( x.size() == 2, "Size must be an array of length 2" );
    return { 2 * x[0], 2 * x[0], x[1] };
}
std::vector<int> CircleFrustum::getLogicalGridSize( const std::vector<double> &res ) const
{
    AMP_INSIST( res.size() == 3u, "Resolution must be an array of length 3" );
    double R = std::max( d_r[0], d_r[1] );
    return { (int) ( R / res[0] ), (int) ( R / res[1] ), (int) ( d_h / res[2] ) };
}


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


/********************************************************
 * Compare the geometry                                  *
 ********************************************************/
bool CircleFrustum::operator==( const Geometry &rhs ) const
{
    auto geom = dynamic_cast<const CircleFrustum *>( &rhs );
    if ( !geom )
        return false;
    return d_dir == geom->d_dir && d_r == geom->d_r && d_h == geom->d_h &&
           d_offset == geom->d_offset && d_C == geom->d_C && d_theta == geom->d_theta;
}


/****************************************************************
 * Write/Read restart data                                       *
 ****************************************************************/
void CircleFrustum::writeRestart( int64_t ) const { AMP_ERROR( "Not finished" ); }
CircleFrustum::CircleFrustum( int64_t ) { AMP_ERROR( "Not finished" ); }


} // namespace AMP::Geometry
