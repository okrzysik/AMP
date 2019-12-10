#include "AMP/ampmesh/shapes/Cylinder.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
Cylinder::Cylinder( AMP::shared_ptr<AMP::Database> db )
{
    d_physicalDim = 3;
    d_logicalDim  = 3;
    d_offset[0]   = 0;
    d_offset[1]   = 0;
    d_offset[2]   = 0;
    auto range    = db->getVector<double>( "Range" );
    AMP_INSIST( range.size() == 3u, "Range must be an array of length 3" );
    d_r     = range[0];
    d_z_min = range[1];
    d_z_max = range[2];
}
Cylinder::Cylinder( double r, double z_min, double z_max )
    : Geometry(), d_r( r ), d_z_min( z_min ), d_z_max( z_max )
{
    d_physicalDim = 3;
    d_logicalDim  = 3;
    d_offset[0]   = 0;
    d_offset[1]   = 0;
    d_offset[2]   = 0;
}


/********************************************************
 * Compute the distance to the object                    *
 * http://mathworld.wolfram.com/Circle-LineIntersection.html
 ********************************************************/
double Cylinder::distance( const Point &pos, const Point &ang ) const
{
    // Get the current point in the reference frame of the cylinder
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2] - 0.5 * ( d_z_min + d_z_max );
    // Compute the distance to the cylinder
    double d = GeometryHelpers::distanceToCylinder( d_r, d_z_max - d_z_min, { x, y, z }, ang );
    return d;
}


/********************************************************
 * Check if the point is inside the geometry             *
 ********************************************************/
bool Cylinder::inside( const Point &pos ) const
{
    double x  = pos.x() - d_offset[0];
    double y  = pos.y() - d_offset[1];
    double z  = pos.z() - d_offset[2];
    double t1 = 1e-12 * d_r * d_r;
    double t2 = 1e-12 * std::max( fabs( d_z_min ), fabs( d_z_max ) );
    double r2 = x * x + y * y;
    bool in_r = r2 <= d_r * d_r + t1;
    bool in_z = z >= d_z_min - t2 && z <= d_z_max + t2;
    return in_r && in_z;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int Cylinder::surface( const Point &pos ) const
{
    double x  = pos.x() - d_offset[0];
    double y  = pos.y() - d_offset[1];
    double z  = pos.z() - d_offset[2];
    double r  = sqrt( x * x + y * y );
    double d1 = std::abs( r - d_r );
    double d2 = std::abs( z - d_z_min );
    double d3 = std::abs( r - d_z_max );
    if ( d1 < std::min( d2, d3 ) )
        return 0;
    else if ( d2 < d3 )
        return 1;
    else
        return 2;
}
Point Cylinder::surfaceNorm( const Point &pos ) const
{
    int s = surface( pos );
    if ( s == 1 ) {
        // -z surface
        return { 0, 0, -1 };
    } else if ( s == 2 ) {
        // -z surface
        return { 0, 0, -1 };
    } else {
        // r
        double x = pos.x() - d_offset[0];
        double y = pos.y() - d_offset[1];
        double n = sqrt( x * x + y * y );
        return { x / n, y / n, 0 };
    }
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point Cylinder::physical( const Point &pos ) const
{
    auto tmp = GeometryHelpers::map_logical_circle( d_r, 2, pos[0], pos[1] );
    double x = tmp.first + d_offset[0];
    double y = tmp.second + d_offset[1];
    double z = d_z_min + pos[2] * ( d_z_max - d_z_min ) + d_offset[2];
    return { x, y, z };
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point Cylinder::logical( const Point &pos ) const
{
    auto tmp =
        GeometryHelpers::map_circle_logical( d_r, 2, pos[0] - d_offset[0], pos[1] - d_offset[1] );
    double z = ( pos[2] - d_z_min - d_offset[2] ) / ( d_z_max - d_z_min );
    return Point( tmp.first, tmp.second, z );
}


/********************************************************
 * Return the centroid and bounding box                  *
 ********************************************************/
Point Cylinder::centroid() const
{
    return { d_offset[0], d_offset[1], d_offset[2] + 0.5 * ( d_z_max + d_z_min ) };
}
std::pair<Point, Point> Cylinder::box() const
{
    Point lb = { d_offset[0] - d_r, d_offset[1] - d_r, d_offset[2] + d_z_min };
    Point ub = { d_offset[0] + d_r, d_offset[1] + d_r, d_offset[2] + d_z_max };
    return { lb, ub };
}


/********************************************************
 * Return the logical grid                               *
 ********************************************************/
std::vector<int> Cylinder::getLogicalGridSize( const std::vector<int> &x ) const
{
    AMP_INSIST( x.size() == 2u, "Size must be an array of length 2" );
    return { x[1], x[1] / 2, x[0] };
}
std::vector<bool> Cylinder::getPeriodicDim() const { return { false, false, false }; }
std::vector<int> Cylinder::getLogicalSurfaceIds() const { return { 4, 4, 4, 4, 2, 1 }; }


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void Cylinder::displaceMesh( const double *x )
{
    d_offset[0] += x[0];
    d_offset[1] += x[1];
    d_offset[2] += x[2];
}


/********************************************************
 * Clone the object                                      *
 ********************************************************/
AMP::shared_ptr<AMP::Geometry::Geometry> Cylinder::clone() const
{
    return AMP::make_shared<Cylinder>( *this );
}


} // namespace Geometry
} // namespace AMP
