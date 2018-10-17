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
    // Initialize the frustrum
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
    // Compute the corners for -/+ x face
    d_face[0][0] = physical( { 0, 0, 0 } );
    d_face[0][1] = physical( { 0, 1, 0 } );
    d_face[0][2] = physical( { 0, 0, 1 } );
    d_face[0][3] = physical( { 0, 1, 1 } );
    d_face[1][0] = physical( { 1, 0, 0 } );
    d_face[1][1] = physical( { 1, 1, 0 } );
    d_face[1][2] = physical( { 1, 0, 1 } );
    d_face[1][3] = physical( { 1, 1, 1 } );
    // Compute the corners for -/+ y face
    d_face[2][0] = physical( { 0, 0, 0 } );
    d_face[2][1] = physical( { 1, 0, 0 } );
    d_face[2][2] = physical( { 0, 0, 1 } );
    d_face[2][3] = physical( { 1, 0, 1 } );
    d_face[3][0] = physical( { 0, 1, 0 } );
    d_face[3][1] = physical( { 1, 1, 0 } );
    d_face[3][2] = physical( { 0, 1, 1 } );
    d_face[3][3] = physical( { 1, 1, 1 } );
    // Compute the corners for -/+ z face
    d_face[4][0] = physical( { 0, 0, 0 } );
    d_face[4][1] = physical( { 1, 0, 0 } );
    d_face[4][2] = physical( { 0, 1, 0 } );
    d_face[4][3] = physical( { 1, 1, 0 } );
    d_face[5][0] = physical( { 0, 0, 1 } );
    d_face[5][1] = physical( { 1, 0, 1 } );
    d_face[5][2] = physical( { 0, 1, 1 } );
    d_face[5][3] = physical( { 1, 1, 1 } );
    // Compute the face normals
    for ( int i = 0; i < 6; i++ )
        d_normal[i] = normal( d_face[i][0], d_face[i][1], d_face[i][2] );
    if ( d_normal[0].x() > 0 )
        d_normal[0] = -d_normal[0];
    if ( d_normal[1].x() < 0 )
        d_normal[1] = -d_normal[1];
    if ( d_normal[2].y() > 0 )
        d_normal[2] = -d_normal[2];
    if ( d_normal[3].y() > 0 )
        d_normal[3] = -d_normal[3];
    if ( d_normal[4].z() > 0 )
        d_normal[4] = -d_normal[4];
    if ( d_normal[5].z() < 0 )
        d_normal[5] = -d_normal[5];
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double SquareFrustum::distance( const Point &pos, const Point &ang ) const
{
    double d = std::numeric_limits<double>::infinity();
    bool in  = inside( pos );
    // Loop over each surface keeping the closest surface
    for ( int i = 0; i < 6; i++ ) {
        // Get the distance from the ray to the plane
        double t = intersectPlane( d_normal[i], d_face[i][0], pos, ang );
        if ( t >= d )
            continue;
        // Check if the point lies within the surface
        auto p = SquareFrustum::logical( pos + t * ang );
        if ( fabs( p[i / 2] ) < 1e-12 || fabs( p[i / 2] - 1.0 ) < 1e-12 ) {
            // The ray intersects the face
            if ( in )
                d = std::min( d, -t );
            else
                d = std::min( d, t );
        }
    }
    return d;
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
    int s     = -1;
    double d2 = std::numeric_limits<double>::infinity();
    for ( int i = 0; i < 6; i++ ) {
        // Get the distance from the ray to the plane
        double t = dot( d_face[i][0] - pos, d_normal[i] );
        // Get the nearest point on the plane
        auto p0 = pos + t * d_normal[i];
        // Restrict the point to the bounds of the surface
        auto p1 = SquareFrustum::logical( p0 );
        p1.x()  = std::min( std::max( p1.x(), 0.0 ), 1.0 );
        p1.y()  = std::min( std::max( p1.y(), 0.0 ), 1.0 );
        p1.z()  = std::min( std::max( p1.z(), 0.0 ), 1.0 );
        auto p  = physical( p1 );
        // Calculate the distance to see if it is closer to the surface
        auto t1   = pos - p;
        double t2 = dot( t1, t1 );
        if ( t2 < d2 ) {
            s  = i;
            d2 = t2;
        }
    }
    return s;
}
Point SquareFrustum::surfaceNorm( const Point &pos ) const
{
    int s = surface( pos );
    return d_normal[s];
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
