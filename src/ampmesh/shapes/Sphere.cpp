#include "AMP/ampmesh/shapes/Sphere.h"
#include "AMP/ampmesh/structured/BoxMeshHelpers.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
Sphere::Sphere( double r ) : d_r( r ) { d_offset.fill( 0 ); }


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double Sphere::distance( const Point<double> &pos, const Point<double> &ang ) const
{
    NULL_USE( pos );
    NULL_USE( ang );
    AMP_ERROR( "Not finished" );
    return 0;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool Sphere::inside( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return false;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int Sphere::surface( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return 0;
}
Point<double> Sphere::surfaceNorm( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point<double>();
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point<double> Sphere::physical( const Point<double> &pos ) const
{
    auto point = AMP::Mesh::BoxMeshHelpers::map_logical_sphere( d_r, pos.x, pos.y, pos.z );
    point.x += d_offset[0];
    point.y += d_offset[1];
    point.z += d_offset[2];
    return point;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point<double> Sphere::logical( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point<double>();
}


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void Sphere::displaceMesh( const double *x )
{
    d_offset[0] += x[0];
    d_offset[1] += x[1];
    d_offset[2] += x[2];
}


} // namespace Geometry
} // namespace AMP
