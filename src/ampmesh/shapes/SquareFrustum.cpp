#include "AMP/ampmesh/shapes/SquareFrustum.h"
#include "AMP/utils/Utilities.h"


#include <vector>


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructors                                          *
 ********************************************************/
SquareFrustum::SquareFrustum( const std::vector<double> &range, int dir, double height )
    : d_dir( dir )
{
    AMP_INSIST( range.size() == 6, "Invalid size for range" );
    AMP_INSIST( dir >= 0 && dir < 6, "Invalid value for dir" );
    for ( size_t i = 0; i < 6; i++ )
        d_range[i] = range[i];
    d_pyramid_size[0] = d_range[1] - d_range[0];
    d_pyramid_size[1] = d_range[3] - d_range[2];
    d_pyramid_size[2] = d_range[5] - d_range[4];
    AMP_INSIST( d_pyramid_size[dir / 2] < height, "Invalid value for height" );
    d_scale_height          = height / d_pyramid_size[dir / 2];
    d_pyramid_size[dir / 2] = height;
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double SquareFrustum::distance( const Point &pos, const Point &ang ) const
{
    NULL_USE( pos );
    NULL_USE( ang );
    AMP_ERROR( "Not finished" );
    return 0;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool SquareFrustum::inside( const Point &pos ) const
{
    // Compute the logical coordinates
    auto p = SquareFrustum::logical( pos );
    return p.x() >= 0 && p.x() <= 1 && p.y() >= 0 && p.y() <= 1 && p.z() >= 0 && p.z() <= 1;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int SquareFrustum::surface( const Point &pos ) const
{
    // Compute the logical coordinates
    // auto p = SquareFrustum::logical( pos );
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return 0;
}
Point SquareFrustum::surfaceNorm( const Point &pos ) const
{
    // Get the surface id
    int s = surface( pos );
    // Set the normal
    Point norm( 0, 0, 0 );
    NULL_USE( s );
    AMP_ERROR( "Not finished" );
    return norm;
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point SquareFrustum::physical( const Point &pos ) const
{
    Point p = pos;
    // Get the point in [0,1,0,1,0,1]
    uint8_t dir2 = d_dir / 2;
    p[d_dir / 2] /= d_scale_height;
    if ( dir2 == 0 ) {
        p = { p.x(), 0.5 + ( p.x() - 1 ) * ( 0.5 - p.y() ), 0.5 + ( p.x() - 1 ) * ( 0.5 - p.z() ) };
    } else if ( dir2 == 1 ) {
        p = { 0.5 + ( p.y() - 1 ) * ( 0.5 - p.x() ), p.y(), 0.5 + ( p.y() - 1 ) * ( 0.5 - p.z() ) };
    } else if ( dir2 == 2 ) {
        p = { 0.5 + ( p.z() - 1 ) * ( 0.5 - p.x() ), 0.5 + ( p.z() - 1 ) * ( 0.5 - p.y() ), p.z() };
    }
    p[dir2] *= d_scale_height;
    if ( d_dir % 2 == 1 )
        p = 1 - p;
    // Get the final coordinates
    p.x() = d_range[0] + p.x() * ( d_range[1] - d_range[0] );
    p.y() = d_range[2] + p.y() * ( d_range[3] - d_range[2] );
    p.z() = d_range[4] + p.z() * ( d_range[5] - d_range[4] );
    return p;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point SquareFrustum::logical( const Point &pos ) const
{
    Point p = pos;
    // Get the point in [0,1,0,1,0,1]
    p.x() = ( p.x() - d_range[0] ) / ( d_range[1] - d_range[0] );
    p.y() = ( p.y() - d_range[2] ) / ( d_range[3] - d_range[2] );
    p.z() = ( p.z() - d_range[4] ) / ( d_range[5] - d_range[4] );
    // Get the final coordinates
    uint8_t dir2 = d_dir / 2;
    if ( d_dir % 2 == 1 )
        p = 1 - p;
    p[dir2] /= d_scale_height;
    if ( dir2 == 0 ) {
        p = { p.x(), 0.5 - ( p.y() - 0.5 ) / ( p.x() - 1 ), 0.5 - ( p.z() - 0.5 ) / ( p.x() - 1 ) };
    } else if ( dir2 == 1 ) {
        p = { 0.5 - ( p.x() - 0.5 ) / ( p.y() - 1 ), p.y(), 0.5 - ( p.z() - 0.5 ) / ( p.y() - 1 ) };
    } else if ( dir2 == 2 ) {
        p = { 0.5 - ( p.x() - 0.5 ) / ( p.z() - 1 ), 0.5 - ( p.y() - 0.5 ) / ( p.z() - 1 ), p.z() };
    }
    p[dir2] *= d_scale_height;
    return p;
}


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void SquareFrustum::displaceMesh( const double *x )
{
    d_range[0] += x[0];
    d_range[1] += x[0];
    d_range[2] += x[1];
    d_range[3] += x[1];
    d_range[4] += x[2];
    d_range[5] += x[2];
}


} // namespace Geometry
} // namespace AMP
