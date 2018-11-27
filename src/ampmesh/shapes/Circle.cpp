#include "AMP/ampmesh/shapes/Circle.h"
#include "AMP/ampmesh/structured/BoxMeshHelpers.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
Circle::Circle( double R ) : d_R( R )
{
    d_offset[0] = 0;
    d_offset[1] = 0;
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double Circle::distance( const Point &pos, const Point &ang ) const
{
    double x  = pos[0] - d_offset[0];
    double y  = pos[1] - d_offset[1];
    double R2 = x * x + y * y;
    if ( R2 < d_R * d_R )
        return sqrt( R2 );
    NULL_USE( ang );
    AMP_ERROR( "Not finished" );
    return 0;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool Circle::inside( const Point &pos ) const
{
    double x = pos[0] - d_offset[0];
    double y = pos[1] - d_offset[1];
    return x * x + y * y <= d_R * d_R;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int Circle::surface( const Point &pos ) const
{
    NULL_USE( pos );
    return 0;
}
Point Circle::surfaceNorm( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point();
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point Circle::physical( const Point &pos ) const
{
    auto tmp = AMP::Mesh::BoxMeshHelpers::map_logical_circle( d_R, 2, pos[0], pos[1] );
    double x = tmp.first + d_offset[0];
    double y = tmp.second + d_offset[1];
    return { x, y };
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point Circle::logical( const Point &pos ) const
{
    auto tmp = AMP::Mesh::BoxMeshHelpers::map_circle_logical( d_R, 2, pos[0], pos[1] );
    return Point( tmp.first, tmp.second );
}


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void Circle::displaceMesh( const double *x )
{
    d_offset[0] += x[0];
    d_offset[1] += x[1];
}


} // namespace Geometry
} // namespace AMP
