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
double Tube::distance( const Point<double> &pos, const Point<double> &ang ) const
{
    NULL_USE( pos );
    NULL_USE( ang );
    AMP_ERROR( "Not finished" );
    return 0;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool Tube::inside( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return false;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int Tube::surface( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return 0;
}
Point<double> Tube::surfaceNorm( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point<double>();
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point<double> Tube::physical( const Point<double> &pos ) const
{
    constexpr double pi = 3.141592653589793116;
    // Compute r, theta
    double r     = d_r_min + pos.x * ( d_r_max - d_r_min );
    double theta = 2.0 * pi * pos.y;
    // Compute the physical coordinate
    Point<double> coord;
    coord.x = r * cos( theta ) + d_offset[0];
    coord.y = r * sin( theta ) + d_offset[1];
    coord.z = d_z_min + pos.z * ( d_z_max - d_z_min ) + d_offset[2];
    return coord;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point<double> Tube::logical( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point<double>();
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
