#include "AMP/ampmesh/shapes/SphereSurface.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
SphereSurface::SphereSurface( std::shared_ptr<AMP::Database> db )
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
SphereSurface::SphereSurface( double r ) : LogicalGeometry(), d_r( r )
{
    d_physicalDim = 3;
    d_logicalDim  = 2;
    d_offset[0]   = 0;
    d_offset[1]   = 0;
    d_offset[2]   = 0;
}


/********************************************************
 * Compute the nearest point on the surface              *
 ********************************************************/
Point SphereSurface::nearest( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return {};
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double SphereSurface::distance( const Point &pos, const Point &ang ) const
{
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2];
    double d = GeometryHelpers::distanceToSphere( d_r, { x, y, z }, ang );
    return std::abs( d );
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool SphereSurface::inside( const Point &pos ) const
{
    double x  = pos.x() - d_offset[0];
    double y  = pos.y() - d_offset[1];
    double z  = pos.z() - d_offset[2];
    double r2 = x * x + y * y + z * z;
    return fabs( r2 - d_r * d_r ) <= 1e-12 * d_r * d_r;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
Point SphereSurface::surfaceNorm( const Point &pos ) const
{
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2];
    double r = sqrt( x * x + y * y + z * z );
    if ( r < d_r )
        return { -x / r, -y / r, -z / r };
    return { x / r, y / r, z / r };
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point SphereSurface::physical( const Point &pos ) const
{
    auto point = GeometryHelpers::map_logical_sphere_surface( d_r, pos[0], pos[1] );
    point[0] += d_offset[0];
    point[1] += d_offset[1];
    point[2] += d_offset[2];
    return point;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point SphereSurface::logical( const Point &pos ) const
{
    double x0 = pos[0] - d_offset[0];
    double y0 = pos[1] - d_offset[1];
    double z0 = pos[2] - d_offset[2];
    auto tmp  = GeometryHelpers::map_sphere_surface_logical( d_r, x0, y0, z0 );
    return Point( tmp.first, tmp.second );
}


/********************************************************
 * Return the centroid and bounding box                  *
 ********************************************************/
Point SphereSurface::centroid() const { return { d_offset[0], d_offset[1], d_offset[2] }; }
std::pair<Point, Point> SphereSurface::box() const
{
    Point lb = { d_offset[0] - d_r, d_offset[1] - d_r, d_offset[2] - d_r };
    Point ub = { d_offset[0] + d_r, d_offset[1] + d_r, d_offset[2] + d_r };
    return { lb, ub };
}


/********************************************************
 * Return the volume                                     *
 ********************************************************/
double SphereSurface::volume() const
{
    constexpr double pi = 3.141592653589793;
    return 4 * pi * d_r * d_r;
}


/********************************************************
 * Return the logical grid                               *
 ********************************************************/
std::vector<int> SphereSurface::getLogicalGridSize( const std::vector<int> &x ) const
{
    AMP_INSIST( x.size() == 1u, "Size must be an array of length 1" );
    return { x[0], x[0] / 2, 1 };
}
std::vector<int> SphereSurface::getLogicalGridSize( const std::vector<double> &res ) const
{
    AMP_INSIST( res.size() == 3u, "Resolution must be an array of length 3" );
    AMP_ERROR( "Not finished" );
    return {};
}
std::vector<bool> SphereSurface::getPeriodicDim() const { return { true, false, false }; }
std::vector<int> SphereSurface::getLogicalSurfaceIds() const { return { -1, -1, -1, -1, -1, -1 }; }


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void SphereSurface::displace( const double *x )
{
    d_offset[0] += x[0];
    d_offset[1] += x[1];
    d_offset[2] += x[2];
}


/********************************************************
 * Clone the object                                      *
 ********************************************************/
std::unique_ptr<AMP::Geometry::Geometry> SphereSurface::clone() const
{
    return std::make_unique<SphereSurface>( *this );
}


/********************************************************
 * Compare the geometry                                  *
 ********************************************************/
bool SphereSurface::operator==( const Geometry &rhs ) const
{
    auto geom = dynamic_cast<const SphereSurface *>( &rhs );
    if ( !geom )
        return false;
    return d_r == geom->d_r && d_offset == geom->d_offset;
}


} // namespace Geometry
} // namespace AMP
