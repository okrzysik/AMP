#include "AMP/ampmesh/shapes/Sphere.h"
#include "AMP/ampmesh/structured/BoxMeshHelpers.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
Sphere::Sphere( double r ) : d_r( r )
{
    d_offset[0] = 0;
    d_offset[1] = 0;
    d_offset[2] = 0;
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double Sphere::distance( const Point &pos, const Point &ang ) const
{
    NULL_USE( pos );
    NULL_USE( ang );
    AMP_ERROR( "Not finished" );
    return 0;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool Sphere::inside( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return false;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int Sphere::surface( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return 0;
}
Point Sphere::surfaceNorm( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point();
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point Sphere::physical( const Point &pos ) const
{
    auto point = AMP::Mesh::BoxMeshHelpers::map_logical_sphere( d_r, pos[0], pos[1], pos[2] );
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
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point();
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
