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
 * Check if the point is inside the geometry             *
 ********************************************************/
bool Tube::inside( const Point &pos ) const
{
    double x  = pos.x() - d_offset[0];
    double y  = pos.y() - d_offset[1];
    double z  = pos.z() - d_offset[2];
    double r2 = x * x + y * y;
    double t1 = 1e-12 * std::max( d_r_min * d_r_min, d_r_max * d_r_max );
    double t2 = 1e-12 * std::max( fabs( d_z_min ), fabs( d_z_max ) );
    bool in_r = r2 >= d_r_min * d_r_min - t1 && r2 <= d_r_max * d_r_max + t1;
    bool in_z = z >= d_z_min - t2 && z <= d_z_max + t2;
    return in_r && in_z;
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
    double theta = 2.0 * pi * ( pos[1] - 0.5 );
    // Compute the physical coordinate
    double x = r * cos( theta ) + d_offset[0];
    double y = r * sin( theta ) + d_offset[1];
    double z = d_z_min + pos[2] * ( d_z_max - d_z_min ) + d_offset[2];
    return { x, y, z };
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point Tube::logical( const Point &pos ) const
{
    constexpr double pi = 3.141592653589793116;
    // Compute r, theta
    double r     = sqrt( ( pos[0] - d_offset[0] ) * ( pos[0] - d_offset[0] ) +
                     ( pos[1] - d_offset[1] ) * ( pos[1] - d_offset[1] ) );
    double theta = acos( ( pos[0] - d_offset[0] ) / r );
    if ( asin( ( pos[1] - d_offset[1] ) / r ) < 0 )
        theta = -theta;
    // Compute the logical coordinate
    double x = ( r - d_r_min ) / ( d_r_max - d_r_min );
    double y = 0.5 + theta / ( 2.0 * pi );
    double z = ( pos[2] - d_z_min - d_offset[2] ) / ( d_z_max - d_z_min );
    return { x, y, z };
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
