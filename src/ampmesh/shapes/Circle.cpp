#include "AMP/ampmesh/shapes/Circle.h"
#include "AMP/ampmesh/structured/BoxMeshHelpers.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
Circle::Circle( double R ) : d_R( R ) { d_offset.fill( 0 ); }


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double Circle::distance( const Point<double> &pos, const Point<double> &ang ) const
{
    double x  = pos.x - d_offset[0];
    double y  = pos.y - d_offset[1];
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
bool Circle::inside( const Point<double> &pos ) const
{
    double x = pos.x - d_offset[0];
    double y = pos.y - d_offset[1];
    return x * x + y * y <= d_R * d_R;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int Circle::surface( const Point<double> &pos ) const
{
    NULL_USE( pos );
    return 0;
}
Point<double> Circle::surfaceNorm( const Point<double> &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return Point<double>();
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point<double> Circle::physical( const Point<double> &pos ) const
{
    auto tmp = AMP::Mesh::BoxMeshHelpers::map_logical_circle( d_R, 2, pos.x, pos.y );
    Point<double> point;
    point.x = tmp.first + d_offset[0];
    point.y = tmp.second + d_offset[1];
    return point;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point<double> Circle::logical( const Point<double> &pos ) const
{
    auto tmp = AMP::Mesh::BoxMeshHelpers::map_circle_logical( d_R, 2, pos.x, pos.y );
    return Point<double>( tmp.first, tmp.second );
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
