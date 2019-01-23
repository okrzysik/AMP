#include "AMP/ampmesh/shapes/Shell.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
Shell::Shell( double r_min, double r_max ) : Geometry(), d_r_min( r_min ), d_r_max( r_max )
{
    d_physicalDim = 3;
    d_logicalDim  = 2;
    d_offset[0]   = 0;
    d_offset[1]   = 0;
    d_offset[2]   = 0;
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
 * Check if the ray is inside the geometry               *
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
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return 0;
}
Point Shell::surfaceNorm( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point();
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
 * Displace the mesh                                     *
 ********************************************************/
void Shell::displaceMesh( const double *x )
{
    d_offset[0] += x[0];
    d_offset[1] += x[1];
    d_offset[2] += x[2];
}


} // namespace Geometry
} // namespace AMP
