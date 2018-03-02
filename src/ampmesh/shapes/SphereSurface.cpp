#include "AMP/ampmesh/shapes/SphereSurface.h"
#include "AMP/ampmesh/structured/BoxMeshHelpers.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
SphereSurface::SphereSurface( double r ) : d_r( r ) { d_offset.fill( 0 ); }


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double SphereSurface::distance( const Point<double> &pos, const Point<double> &ang ) const
{
    NULL_USE( pos );
    NULL_USE( ang );
    AMP_ERROR( "Not finished" );
    return 0;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool SphereSurface::inside( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return false;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int SphereSurface::surface( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return 0;
}
Point<double> SphereSurface::surfaceNorm( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point<double>();
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point<double> SphereSurface::physical( const Point<double> &pos ) const
{
    auto point = AMP::Mesh::BoxMeshHelpers::map_logical_sphere_surface( d_r, pos.x, pos.y );
    point.x += d_offset[0];
    point.y += d_offset[1];
    point.z += d_offset[2];
    return point;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point<double> SphereSurface::logical( const Point<double> &pos ) const
{
    double x0 = pos.x - d_offset[0];
    double y0 = pos.y - d_offset[1];
    double z0 = pos.z - d_offset[2];
    double r  = sqrt( x0 * x0 + y0 * y0 + z0 * z0 );
    auto tmp  = AMP::Mesh::BoxMeshHelpers::map_sphere_surface_logical( d_r, x0, y0, z0 );
    return Point<double>( tmp.first, tmp.second, r / d_r - 1 );
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
