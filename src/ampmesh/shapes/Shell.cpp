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
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return false;
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
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point();
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
