#include "AMP/ampmesh/shapes/Sphere.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
Sphere::Sphere( std::shared_ptr<AMP::Database> db )
{
    d_physicalDim = 3;
    d_logicalDim  = 2;
    d_offset[0]   = 0;
    d_offset[1]   = 0;
    d_offset[2]   = 0;
    auto range    = db->getVector<double>( "Range" );
    AMP_INSIST( range.size() == 1u, "Range must be an array of length 1" );
    d_r = range[0];
}
Sphere::Sphere( double r ) : LogicalGeometry(), d_r( r )
{
    d_physicalDim = 3;
    d_logicalDim  = 3;
    d_offset[0]   = 0;
    d_offset[1]   = 0;
    d_offset[2]   = 0;
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double Sphere::distance( const Point &pos, const Point &ang ) const
{
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2];
    double d = GeometryHelpers::distanceToSphere( d_r, { x, y, z }, ang );
    return d;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool Sphere::inside( const Point &pos ) const
{
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2];
    return x * x + y * y + z * z <= ( 1.0 + 1e-12 ) * d_r * d_r;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
Point Sphere::surfaceNorm( const Point &pos ) const
{
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2];
    double r = sqrt( x * x + y * y + z * z );
    return { x / r, y / r, z / r };
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point Sphere::physical( const Point &pos ) const
{
    auto point = GeometryHelpers::map_logical_sphere( d_r, pos.x(), pos.y(), pos.z() );
    point[0] += d_offset[0];
    point[1] += d_offset[1];
    point[2] += d_offset[2];
    return point;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point Sphere::logical( const Point &pos ) const
{
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2];
    return GeometryHelpers::map_sphere_logical( d_r, x, y, z );
}


/********************************************************
 * Return the centroid and bounding box                  *
 ********************************************************/
Point Sphere::centroid() const { return { d_offset[0], d_offset[1], d_offset[2] }; }
std::pair<Point, Point> Sphere::box() const
{
    Point lb = { d_offset[0] - d_r, d_offset[1] - d_r, d_offset[2] - d_r };
    Point ub = { d_offset[0] + d_r, d_offset[1] + d_r, d_offset[2] + d_r };
    return { lb, ub };
}


/********************************************************
 * Return the logical grid                               *
 ********************************************************/
std::vector<int> Sphere::getLogicalGridSize( const std::vector<int> &x ) const
{
    AMP_INSIST( x.size() == 1u, "Size must be an array of length 1" );
    return { 2 * x[0], 2 * x[0], 2 * x[0] };
}
std::vector<bool> Sphere::getPeriodicDim() const { return { false, false, false }; }
std::vector<int> Sphere::getLogicalSurfaceIds() const { return { 4, 4, 4, 4, 2, 1 }; }


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void Sphere::displace( const double *x )
{
    d_offset[0] += x[0];
    d_offset[1] += x[1];
    d_offset[2] += x[2];
}


/********************************************************
 * Clone the object                                      *
 ********************************************************/
std::unique_ptr<AMP::Geometry::Geometry> Sphere::clone() const
{
    return std::make_unique<Sphere>( *this );
}


} // namespace Geometry
} // namespace AMP
