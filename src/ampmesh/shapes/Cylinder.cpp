#include "AMP/ampmesh/shapes/Cylinder.h"
#include "AMP/ampmesh/structured/BoxMeshHelpers.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
Cylinder::Cylinder( double r, double z_min, double z_max )
    : d_r( r ), d_z_min( z_min ), d_z_max( z_max )
{
    d_offset[0] = 0;
    d_offset[1] = 0;
    d_offset[2] = 0;
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double Cylinder::distance( const Point &pos, const Point &ang ) const
{
    NULL_USE( pos );
    NULL_USE( ang );
    AMP_ERROR( "Not finished" );
    return 0;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool Cylinder::inside( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return false;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int Cylinder::surface( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return 0;
}
Point Cylinder::surfaceNorm( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point();
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point Cylinder::physical( const Point &pos ) const
{
    auto tmp = AMP::Mesh::BoxMeshHelpers::map_logical_circle( d_r, 2, pos[0], pos[1] );
    Point point;
    point[0] = tmp.first + d_offset[0];
    point[1] = tmp.second + d_offset[1];
    point[2] = d_z_min + pos[2] * ( d_z_max - d_z_min ) + d_offset[2];
    return point;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point Cylinder::logical( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point();
}


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void Cylinder::displaceMesh( const double *x )
{
    d_offset[0] += x[0];
    d_offset[1] += x[1];
    d_offset[2] += x[2];
}


} // namespace Geometry
} // namespace AMP
