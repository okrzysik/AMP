#include "AMP/ampmesh/shapes/Shell.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
Shell::Shell( std::shared_ptr<const AMP::Database> db )
{
    d_physicalDim = 3;
    d_logicalDim  = 3;
    d_offset[0]   = 0;
    d_offset[1]   = 0;
    d_offset[2]   = 0;
    auto range    = db->getVector<double>( "Range" );
    AMP_INSIST( range.size() == 2u, "Range must be an array of length 2" );
    d_r_min = range[0];
    d_r_max = range[1];
}
Shell::Shell( double r_min, double r_max ) : LogicalGeometry(), d_r_min( r_min ), d_r_max( r_max )
{
    d_physicalDim = 3;
    d_logicalDim  = 3;
    d_offset[0]   = 0;
    d_offset[1]   = 0;
    d_offset[2]   = 0;
}


/********************************************************
 * Compute the nearest point on the surface              *
 ********************************************************/
Point Shell::nearest( const Point &pos ) const
{
    // Get the current point in the reference frame of the circle
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2];
    // Calculate the nearest point
    double r = sqrt( x * x + y * y + z * z );
    if ( r == 0 ) {
        x = d_r_min;
    } else if ( r < d_r_min ) {
        x *= d_r_min / r;
        y *= d_r_min / r;
        z *= d_r_min / r;
    } else if ( r > d_r_max ) {
        x *= d_r_max / r;
        y *= d_r_max / r;
        z *= d_r_max / r;
    }
    return { x + d_offset[0], y + d_offset[1], z + d_offset[2] };
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double Shell::distance( const Point &pos, const Point &ang ) const
{
    double x    = pos.x() - d_offset[0];
    double y    = pos.y() - d_offset[1];
    double z    = pos.z() - d_offset[2];
    double d1   = GeometryHelpers::distanceToSphere( d_r_min, { x, y, z }, ang );
    double d2   = GeometryHelpers::distanceToSphere( d_r_max, { x, y, z }, ang );
    double d    = std::min( std::abs( d1 ), std::abs( d2 ) );
    double r2   = x * x + y * y + z * z;
    bool inside = r2 >= d_r_min * d_r_min && r2 <= d_r_max * d_r_max;
    return ( inside ? -1 : 1 ) * d;
}


/********************************************************
 * Check if the point is inside the geometry             *
 ********************************************************/
bool Shell::inside( const Point &pos ) const
{
    double x      = pos.x() - d_offset[0];
    double y      = pos.y() - d_offset[1];
    double z      = pos.z() - d_offset[2];
    double r2     = x * x + y * y + z * z;
    double r2_min = ( 1.0 - 1e-12 ) * d_r_min * d_r_min;
    double r2_max = ( 1.0 + 1e-12 ) * d_r_max * d_r_max;
    return r2 >= r2_min && r2 <= r2_max;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int Shell::surface( const Point &pos ) const
{
    double x  = pos.x() - d_offset[0];
    double y  = pos.y() - d_offset[1];
    double z  = pos.z() - d_offset[2];
    double r2 = x * x + y * y + z * z;
    bool test = fabs( r2 - d_r_min * d_r_min ) < fabs( r2 - d_r_max * d_r_max );
    return test ? 0 : 1;
}
Point Shell::surfaceNorm( const Point &pos ) const
{
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2];
    double r = sqrt( x * x + y * y + z * z );
    if ( fabs( r - d_r_min ) < fabs( r - d_r_max ) )
        return { -x / r, -y / r, -z / r };
    return { x / r, y / r, z / r };
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point Shell::physical( const Point &pos ) const
{
    auto point = GeometryHelpers::map_logical_shell( d_r_min, d_r_max, pos[0], pos[1], pos[2] );
    point[0] += d_offset[0];
    point[1] += d_offset[1];
    point[2] += d_offset[2];
    return point;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point Shell::logical( const Point &pos ) const
{
    return GeometryHelpers::map_shell_logical(
        d_r_min, d_r_max, pos[0] - d_offset[0], pos[1] - d_offset[1], pos[2] - d_offset[2] );
}


/********************************************************
 * Return the centroid and bounding box                  *
 ********************************************************/
Point Shell::centroid() const { return { d_offset[0], d_offset[1], d_offset[2] }; }
std::pair<Point, Point> Shell::box() const
{
    Point lb = { d_offset[0] - d_r_max, d_offset[1] - d_r_max, d_offset[2] - d_r_max };
    Point ub = { d_offset[0] + d_r_max, d_offset[1] + d_r_max, d_offset[2] + d_r_max };
    return { lb, ub };
}


/********************************************************
 * Return the volume                                     *
 ********************************************************/
double Shell::volume() const
{
    constexpr double pi = 3.141592653589793;
    double V1           = 4.0 / 3.0 * pi * d_r_max * d_r_max * d_r_max;
    double V2           = 4.0 / 3.0 * pi * d_r_min * d_r_min * d_r_min;
    return V1 - V2;
}


/********************************************************
 * Return the logical grid                               *
 ********************************************************/
std::vector<int> Shell::getLogicalGridSize( const std::vector<int> &x ) const
{
    AMP_INSIST( x.size() == 2u, "Size must be an array of length 2" );
    return { 2 * x[0], 2 * x[0], x[1] };
}
std::vector<int> Shell::getLogicalGridSize( const std::vector<double> &res ) const
{
    AMP_INSIST( res.size() == 3u, "Resolution must be an array of length 3" );
    AMP_ERROR( "Not finished" );
    return {};
}
std::vector<bool> Shell::getPeriodicDim() const { return { true, false, false }; }
std::vector<int> Shell::getLogicalSurfaceIds() const { return { -1, -1, 1, 2, 3, 4 }; }


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void Shell::displace( const double *x )
{
    d_offset[0] += x[0];
    d_offset[1] += x[1];
    d_offset[2] += x[2];
}


/********************************************************
 * Clone the object                                      *
 ********************************************************/
std::unique_ptr<AMP::Geometry::Geometry> Shell::clone() const
{
    return std::make_unique<Shell>( *this );
}


/********************************************************
 * Compare the geometry                                  *
 ********************************************************/
bool Shell::operator==( const Geometry &rhs ) const
{
    auto geom = dynamic_cast<const Shell *>( &rhs );
    if ( !geom )
        return false;
    return d_r_min == geom->d_r_min && d_r_max == geom->d_r_max && d_offset == geom->d_offset;
}


} // namespace Geometry
} // namespace AMP
