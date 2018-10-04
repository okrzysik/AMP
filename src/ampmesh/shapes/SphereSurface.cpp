#include "AMP/ampmesh/shapes/SphereSurface.h"
#include "AMP/ampmesh/structured/BoxMeshHelpers.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
SphereSurface::SphereSurface( double r ) : d_r( r )
{
    d_offset[0] = 0;
    d_offset[1] = 0;
    d_offset[2] = 0;
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double SphereSurface::distance( const Point &pos, const Point &ang ) const
{
    NULL_USE( pos );
    NULL_USE( ang );
    AMP_ERROR( "Not finished" );
    return 0;
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
int SphereSurface::surface( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return 0;
}
Point SphereSurface::surfaceNorm( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point();
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point SphereSurface::physical( const Point &pos ) const
{
    auto point = AMP::Mesh::BoxMeshHelpers::map_logical_sphere_surface( d_r, pos[0], pos[1] );
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
    double r  = sqrt( x0 * x0 + y0 * y0 + z0 * z0 );
    auto tmp  = AMP::Mesh::BoxMeshHelpers::map_sphere_surface_logical( d_r, x0, y0, z0 );
    return Point( tmp.first, tmp.second, r / d_r - 1 );
}


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void SphereSurface::displaceMesh( const double *x )
{
    d_offset[0] += x[0];
    d_offset[1] += x[1];
    d_offset[2] += x[2];
}


} // namespace Geometry
} // namespace AMP
