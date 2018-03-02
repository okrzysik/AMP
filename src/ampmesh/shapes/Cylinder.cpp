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
    d_offset.fill( 0 );
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double Cylinder::distance( const Point<double> &pos, const Point<double> &ang ) const
{
    NULL_USE( pos );
    NULL_USE( ang );
    AMP_ERROR( "Not finished" );
    return 0;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool Cylinder::inside( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return false;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int Cylinder::surface( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return 0;
}
Point<double> Cylinder::surfaceNorm( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point<double>();
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point<double> Cylinder::physical( const Point<double> &pos ) const
{
    auto tmp = AMP::Mesh::BoxMeshHelpers::map_logical_circle( d_r, 2, pos.x, pos.y );
    Point<double> point;
    point.x = tmp.first + d_offset[0];
    point.y = tmp.second + d_offset[1];
    point.z = d_z_min + pos.z * ( d_z_max - d_z_min ) + d_offset[2];
    return point;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point<double> Cylinder::logical( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point<double>();
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
