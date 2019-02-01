#include "AMP/ampmesh/shapes/Tube.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
Tube::Tube( double r_min, double r_max, double z_min, double z_max )
    : Geometry(), d_r_min( r_min ), d_r_max( r_max ), d_z_min( z_min ), d_z_max( z_max )
{
    d_physicalDim = 3;
    d_logicalDim  = 2;
    d_offset[0]   = 0;
    d_offset[1]   = 0;
    d_offset[2]   = 0;
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double Tube::distance( const Point &pos, const Point &ang ) const
{
    // Get the current point in the reference frame of the cylinder
    double h = d_z_max - d_z_min;
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2] - 0.5 * ( d_z_min + d_z_max );
    // Compute the distance to the cylinders
    double d1 = GeometryHelpers::distanceToCylinder( d_r_min, h, { x, y, z }, ang );
    double d2 = GeometryHelpers::distanceToCylinder( d_r_max, h, { x, y, z }, ang );
    double d  = std::min( std::abs( d1 ), std::abs( d2 ) );
    // Compute the point of intersection
    auto p    = Point( { x, y, z } ) + d * ang;
    double r2 = p.x() * p.x() + p.y() * p.y();
    if ( r2 >= d_r_min * d_r_min - 1e-12 ) {
        bool inside = r2 <= d_r_max * d_r_max && std::abs( p.z() ) <= 0.5 * h;
        return ( inside ? -1 : 1 ) * d;
    }
    if ( p.z() * ang.z() > 0 )
        return std::numeric_limits<double>::infinity();
    double d3 = GeometryHelpers::distanceToCylinder( d_r_min, 1e200, p, ang );
    d         = d + std::abs( d3 );
    p         = Point( { x, y, z } ) + d * ang;
    r2        = p.x() * p.x() + p.y() * p.y();
    if ( r2 >= d_r_min * d_r_min - 1e-12 && std::abs( p.z() ) <= 0.5 * h )
        return d;
    return std::numeric_limits<double>::infinity();
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
 * Return the centroid and bounding box                  *
 ********************************************************/
Point Tube::centroid() const
{
    return { d_offset[0], d_offset[1], d_offset[2] + 0.5 * ( d_z_max + d_z_min ) };
}
std::pair<Point, Point> Tube::box() const
{
    Point lb = { d_offset[0] - d_r_max, d_offset[1] - d_r_max, d_offset[2] + d_z_min };
    Point ub = { d_offset[0] + d_r_max, d_offset[1] + d_r_max, d_offset[2] + d_z_max };
    return { lb, ub };
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
