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
 * http://mathworld.wolfram.com/Circle-LineIntersection.html
 ********************************************************/
double Cylinder::distance( const Point &pos, const Point &ang ) const
{
    // Get the current point in the reference frame of the cylinder
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2];
    // Compute the intersection of a line with the circle of the cylinder
    double dx = ang.x();
    double dy = ang.y();
    double dr = sqrt( dx * dx + dy * dy );
    double D  = x * ( y + ang.y() ) - ( x + ang.x() ) * y;
    double t  = d_r * d_r * dr * dr - D * D;
    if ( t < 0 )
        return std::numeric_limits<double>::infinity();
    t         = sqrt( t );
    double s  = dy < 0 ? -1 : 1;
    double x1 = ( D * dy + s * dx * t ) / ( dr * dr );
    double x2 = ( D * dy - s * dx * t ) / ( dr * dr );
    // double y1 = ( -D * dx + abs( dy ) *t ) / ( dr*dr);
    // double y2 = ( -D * dx - abs( dy ) *t ) / ( dr*dr);
    // Compute the distance to the point
    double d1 = ( x1 - x ) / ang.x();
    double d2 = ( x2 - x ) / ang.x();
    if ( d1 < 0 )
        d1 = std::numeric_limits<double>::infinity();
    if ( d2 < 0 )
        d2 = std::numeric_limits<double>::infinity();
    // Check that the z-point is within the cylinder for each point
    double z1 = z + d1 * ang.z();
    double z2 = z + d2 * ang.z();
    if ( z1 < d_z_min || z1 > d_z_max )
        d1 = std::numeric_limits<double>::infinity();
    if ( z2 < d_z_min || z2 > d_z_max )
        d2 = std::numeric_limits<double>::infinity();
    // Return the distance to the closest point
    bool inside = d1 < 1e100 && d2 < 1e100;
    if ( d1 < 1e100 && d2 < 1e100 )
        return std::min( d1, d2 );
    return ( inside ? -1 : 1 ) * std::min( d1, d2 );
    ;
}


/********************************************************
 * Check if the point is inside the geometry             *
 ********************************************************/
bool Cylinder::inside( const Point &pos ) const
{
    double x  = pos.x() - d_offset[0];
    double y  = pos.y() - d_offset[1];
    double z  = pos.z() - d_offset[2];
    double t1 = 1e-12 * d_r * d_r;
    double t2 = 1e-12 * std::max( fabs( d_z_min ), fabs( d_z_max ) );
    double r2 = x * x + y * y;
    bool in_r = r2 <= d_r * d_r + t1;
    bool in_z = z >= d_z_min - t2 && z <= d_z_max + t2;
    return in_r && in_z;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int Cylinder::surface( const Point &pos ) const
{
    double x  = pos.x() - d_offset[0];
    double y  = pos.y() - d_offset[1];
    double z  = pos.z() - d_offset[2];
    double r  = sqrt( x * x + y * y );
    double d1 = std::abs( r - d_r );
    double d2 = std::abs( z - d_z_min );
    double d3 = std::abs( r - d_z_max );
    if ( d1 < std::min( d2, d3 ) )
        return 0;
    else if ( d2 < d3 )
        return 1;
    else
        return 2;
}
Point Cylinder::surfaceNorm( const Point &pos ) const
{
    int s = surface( pos );
    if ( s == 1 ) {
        // -z surface
        return { 0, 0, -1 };
    } else if ( s == 2 ) {
        // -z surface
        return { 0, 0, -1 };
    } else {
        // r
        double x = pos.x() - d_offset[0];
        double y = pos.y() - d_offset[1];
        double n = sqrt( x * x + y * y );
        return { x / n, y / n, 0 };
    }
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point Cylinder::physical( const Point &pos ) const
{
    auto tmp = AMP::Mesh::BoxMeshHelpers::map_logical_circle( d_r, 2, pos[0], pos[1] );
    double x = tmp.first + d_offset[0];
    double y = tmp.second + d_offset[1];
    double z = d_z_min + pos[2] * ( d_z_max - d_z_min ) + d_offset[2];
    return { x, y, z };
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point Cylinder::logical( const Point &pos ) const
{
    auto tmp = AMP::Mesh::BoxMeshHelpers::map_circle_logical(
        d_r, 2, pos[0] - d_offset[0], pos[1] - d_offset[1] );
    double z = ( pos[2] - d_z_min - d_offset[2] ) / ( d_z_max - d_z_min );
    return Point( tmp.first, tmp.second, z );
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
