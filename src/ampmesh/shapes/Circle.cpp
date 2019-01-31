#include "AMP/ampmesh/shapes/Circle.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
Circle::Circle( double R ) : Geometry(), d_R( R )
{
    d_physicalDim = 2;
    d_logicalDim  = 2;
    d_offset[0]   = 0;
    d_offset[1]   = 0;
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double Circle::distance( const Point &pos, const Point &ang ) const
{
    // Get the current point in the reference frame of the cylinder
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = 0;
    // Compute the distance to the cylinder
    double d = GeometryHelpers::distanceToCylinder( d_R, 1e200, { x, y, z }, ang );
    return d;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool Circle::inside( const Point &pos ) const
{
    double x   = pos[0] - d_offset[0];
    double y   = pos[1] - d_offset[1];
    double R21 = x * x + y * y;
    double R22 = d_R * d_R;
    return R21 <= ( 1.0 + 1e-12 ) * R22;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
Point Circle::surfaceNorm( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point();
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point Circle::physical( const Point &pos ) const
{
    auto tmp = GeometryHelpers::map_logical_circle( d_R, 2, pos[0], pos[1] );
    double x = tmp.first + d_offset[0];
    double y = tmp.second + d_offset[1];
    return { x, y };
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point Circle::logical( const Point &pos ) const
{
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    auto tmp = GeometryHelpers::map_circle_logical( d_R, 2, x, y );
    return Point( tmp.first, tmp.second );
}


/********************************************************
 * Return the centroid and bounding box                  *
 ********************************************************/
Point Circle::centroid() const { return { d_offset[0], d_offset[1] }; }
std::pair<Point, Point> Circle::box() const
{
    Point lb = { d_offset[0] - d_R, d_offset[1] - d_R };
    Point ub = { d_offset[0] + d_R, d_offset[1] + d_R };
    return { lb, ub };
}


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void Circle::displaceMesh( const double *x )
{
    d_offset[0] += x[0];
    d_offset[1] += x[1];
}


} // namespace Geometry
} // namespace AMP
