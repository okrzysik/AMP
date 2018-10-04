#include "AMP/ampmesh/shapes/Tube.h"
#include "AMP/ampmesh/structured/BoxMeshHelpers.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
Tube::Tube( double r_min, double r_max, double z_min, double z_max )
    : d_r_min( r_min ), d_r_max( r_max ), d_z_min( z_min ), d_z_max( z_max )
{
    d_offset[0] = 0;
    d_offset[1] = 0;
    d_offset[2] = 0;
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double Tube::distance( const Point &pos, const Point &ang ) const
{
    NULL_USE( pos );
    NULL_USE( ang );
    AMP_ERROR( "Not finished" );
    return 0;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool Tube::inside( const Point &pos ) const
{
    double x  = pos.x() - d_offset[0];
    double y  = pos.y() - d_offset[1];
    double z  = pos.z() - d_offset[2];
    double r2 = x * x + y * y;
    return r2 >= d_r_min * d_r_min && r2 <= d_r_max * d_r_max && z >= d_z_min && z <= d_z_max;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int Tube::surface( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return 0;
}
Point Tube::surfaceNorm( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point();
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point Tube::physical( const Point &pos ) const
{
    constexpr double pi = 3.141592653589793116;
    // Compute r, theta
    double r     = d_r_min + pos[0] * ( d_r_max - d_r_min );
    double theta = 2.0 * pi * pos[1];
    // Compute the physical coordinate
    Point coord;
    coord[0] = r * cos( theta ) + d_offset[0];
    coord[1] = r * sin( theta ) + d_offset[1];
    coord[2] = d_z_min + pos[2] * ( d_z_max - d_z_min ) + d_offset[2];
    return coord;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point Tube::logical( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point();
}


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void Tube::displaceMesh( const double *x )
{
    d_offset[0] += x[0];
    d_offset[1] += x[1];
    d_offset[2] += x[2];
}


} // namespace Geometry
} // namespace AMP
