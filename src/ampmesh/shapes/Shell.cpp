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
    d_offset.fill( 0 );
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double Shell::distance( const Point<double> &pos, const Point<double> &ang ) const
{
    NULL_USE( pos );
    NULL_USE( ang );
    AMP_ERROR( "Not finished" );
    return 0;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool Shell::inside( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return false;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int Shell::surface( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return 0;
}
Point<double> Shell::surfaceNorm( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point<double>();
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point<double> Shell::physical( const Point<double> &pos ) const
{
    auto point =
        AMP::Mesh::BoxMeshHelpers::map_logical_shell( d_r_min, d_r_max, pos.x, pos.y, pos.z );
    point.x += d_offset[0];
    point.y += d_offset[1];
    point.z += d_offset[2];
    return point;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point<double> Shell::logical( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point<double>();
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
