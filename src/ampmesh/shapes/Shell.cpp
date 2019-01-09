#include "AMP/ampmesh/shapes/Shell.h"
#include "AMP/ampmesh/structured/BoxMeshHelpers.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
Shell::Shell( double r_min, double r_max ) : d_r_min( r_min ), d_r_max( r_max )
{
    d_offset[0] = 0;
    d_offset[1] = 0;
    d_offset[2] = 0;
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double Shell::distance( const Point &pos, const Point &ang ) const
{
    NULL_USE( pos );
    NULL_USE( ang );
    AMP_ERROR( "Not finished" );
    return 0;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool Shell::inside( const Point &pos ) const
{
    double x      = pos.x() - d_offset[0];
    double y      = pos.y() - d_offset[1];
    double z      = pos.z() - d_offset[2];
    double r2     = x * x + y * y + z * z;
    double r2_min = ( 1.0 - 1e-12 ) * d_r_min * d_r_min;
    double r2_max = ( 1.0 + 1e-12 ) * d_r_max * d_r_max;
    return r2 >= r2_min && r2 <= r2_max;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int Shell::surface( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return 0;
}
Point Shell::surfaceNorm( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point();
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point Shell::physical( const Point &pos ) const
{
    auto point =
        AMP::Mesh::BoxMeshHelpers::map_logical_shell( d_r_min, d_r_max, pos[0], pos[1], pos[2] );
    point[0] += d_offset[0];
    point[1] += d_offset[1];
    point[2] += d_offset[2];
    return point;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point Shell::logical( const Point &pos ) const
{
    return AMP::Mesh::BoxMeshHelpers::map_shell_logical(
        d_r_min, d_r_max, pos[0] - d_offset[0], pos[1] - d_offset[1], pos[2] - d_offset[2] );
}


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void Shell::displaceMesh( const double *x )
{
    d_offset[0] += x[0];
    d_offset[1] += x[1];
    d_offset[2] += x[2];
}


} // namespace Geometry
} // namespace AMP
