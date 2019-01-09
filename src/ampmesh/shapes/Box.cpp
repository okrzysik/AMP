#include "AMP/ampmesh/shapes/Box.h"
#include "AMP/utils/Utilities.h"


#include <vector>


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructors                                          *
 ********************************************************/
template<std::size_t NDIM>
Box<NDIM>::Box( const std::vector<double> &range )
{
    static_assert( NDIM >= 0 && NDIM <= 3, "Invalid number of dimensions" );
    for ( int i = 0; i < 3; i++ ) {
        d_range[2 * i + 0] = -1e100;
        d_range[2 * i + 1] = 1e100;
    }
    for ( size_t i = 0; i < 2 * NDIM; i++ )
        d_range[i] = range[i];
}
template<std::size_t NDIM>
Grid<NDIM>::Grid( const std::vector<std::vector<double>> &coord )
{
    static_assert( NDIM >= 0 && NDIM <= 3, "Invalid number of dimensions" );
    for ( int i = 0; i < 3; i++ ) {
        d_range[2 * i + 0] = -1e100;
        d_range[2 * i + 1] = 1e100;
    }
    for ( size_t i = 0; i < coord.size(); i++ ) {
        d_coord[i] = coord[i];
        for ( const auto tmp : d_coord[i] ) {
            d_range[2 * i + 0] = std::min( d_range[2 * i + 0], tmp );
            d_range[2 * i + 1] = std::max( d_range[2 * i + 1], tmp );
        }
    }
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
static inline double calcDist( const Point &pos, const Point &ang, const double *range )
{
    constexpr double tol = 1e-5;
    // Compute the distance to each surface
    double d1 = ( range[0] - pos[0] ) / ang[0];
    double d2 = ( range[1] - pos[0] ) / ang[0];
    double d3 = ( range[2] - pos[1] ) / ang[1];
    double d4 = ( range[3] - pos[1] ) / ang[1];
    double d5 = ( range[4] - pos[2] ) / ang[2];
    double d6 = ( range[5] - pos[2] ) / ang[2];
    // Check if point is inside volume (within the tolerance)
    bool inside = ( ang[0] >= 0 ? ( d1 < tol && d2 > tol ) : d2 < tol && d1 > tol ) &&
                  ( ang[1] >= 0 ? ( d3 < tol && d4 > tol ) : d4 < tol && d3 > tol ) &&
                  ( ang[2] >= 0 ? ( d5 < tol && d6 > tol ) : d6 < tol && d5 > tol );
    // Compute the distance
    double dist = std::numeric_limits<double>::quiet_NaN();
    if ( inside ) {
        dist = std::max( d1, d2 );
        dist = std::min( dist, std::max( d3, d4 ) );
        dist = std::min( dist, std::max( d5, d6 ) );
        dist = -dist;
    } else {
        double distx = std::numeric_limits<double>::infinity();
        double disty = std::numeric_limits<double>::infinity();
        double distz = std::numeric_limits<double>::infinity();
        distx        = d1 > 0 ? std::min( distx, d1 ) : distx;
        distx        = d2 > 0 ? std::min( distx, d2 ) : distx;
        disty        = d3 > 0 ? std::min( disty, d3 ) : disty;
        disty        = d4 > 0 ? std::min( disty, d4 ) : disty;
        distz        = d5 > 0 ? std::min( distz, d5 ) : distz;
        distz        = d6 > 0 ? std::min( distz, d6 ) : distz;
        dist         = std::max( distx, std::max( disty, distz ) );
    }
    return dist;
}
template<std::size_t NDIM>
double Box<NDIM>::distance( const Point &pos, const Point &ang ) const
{
    return calcDist( pos, ang, d_range );
}
template<std::size_t NDIM>
double Grid<NDIM>::distance( const Point &pos, const Point &ang ) const
{
    return calcDist( pos, ang, d_range );
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
template<>
bool Box<1>::inside( const Point &pos ) const
{
    return pos[0] >= d_range[0] && pos[0] <= d_range[1];
}
template<>
bool Box<2>::inside( const Point &pos ) const
{
    return pos[0] >= d_range[0] && pos[0] <= d_range[1] && pos[1] >= d_range[2] &&
           pos[1] <= d_range[3];
}
template<>
bool Box<3>::inside( const Point &pos ) const
{
    return pos[0] >= d_range[0] && pos[0] <= d_range[1] && pos[1] >= d_range[2] &&
           pos[1] <= d_range[3] && pos[2] >= d_range[4] && pos[2] <= d_range[5];
}
template<>
bool Grid<1>::inside( const Point &pos ) const
{
    return pos[0] >= d_range[0] && pos[0] <= d_range[1] && pos[1];
}
template<>
bool Grid<2>::inside( const Point &pos ) const
{
    return pos[0] >= d_range[0] && pos[0] <= d_range[1] && pos[1] >= d_range[2] &&
           pos[1] <= d_range[3];
}
template<>
bool Grid<3>::inside( const Point &pos ) const
{
    return pos[0] >= d_range[0] && pos[0] <= d_range[1] && pos[1] >= d_range[2] &&
           pos[1] <= d_range[3] && pos[2] >= d_range[4] && pos[2] <= d_range[5];
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
template<std::size_t NDIM>
int Box<NDIM>::surface( const Point &pos ) const
{
    double d[6] = { fabs( pos[0] - d_range[0] ), fabs( pos[0] - d_range[1] ),
                    fabs( pos[1] - d_range[2] ), fabs( pos[1] - d_range[3] ),
                    fabs( pos[2] - d_range[4] ), fabs( pos[2] - d_range[5] ) };
    double d0   = d[0];
    int s       = 0;
    for ( int i = 1; i < 6; i++ ) {
        if ( d[i] < d0 ) {
            d0 = d[i];
            s  = i;
        }
    }
    return s;
}
template<std::size_t NDIM>
int Grid<NDIM>::surface( const Point &pos ) const
{
    double d[6] = { fabs( pos[0] - d_range[0] ), fabs( pos[0] - d_range[1] ),
                    fabs( pos[1] - d_range[2] ), fabs( pos[1] - d_range[3] ),
                    fabs( pos[2] - d_range[4] ), fabs( pos[2] - d_range[5] ) };
    double d0   = d[0];
    int s       = 0;
    for ( int i = 1; i < 6; i++ ) {
        if ( d[i] < d0 ) {
            d0 = d[i];
            s  = i;
        }
    }
    return s;
}
template<std::size_t NDIM>
Point Box<NDIM>::surfaceNorm( const Point &pos ) const
{
    // Get the surface id
    int s = surface( pos );
    // Set the normal
    Point norm( NDIM );
    norm[s / 2] = s % 2 == 0 ? -1 : 1;
    return norm;
}
template<std::size_t NDIM>
Point Grid<NDIM>::surfaceNorm( const Point &pos ) const
{
    // Get the surface id
    int s = surface( pos );
    // Set the normal
    Point norm( NDIM );
    norm[s / 2] = s % 2 == 0 ? -1 : 1;
    return norm;
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
template<>
Point Box<1>::physical( const Point &pos ) const
{
    return Point( d_range[0] + pos[0] * ( d_range[1] - d_range[0] ) );
}
template<>
Point Box<2>::physical( const Point &pos ) const
{
    return Point( d_range[0] + pos[0] * ( d_range[1] - d_range[0] ),
                  d_range[2] + pos[1] * ( d_range[3] - d_range[2] ) );
}
template<>
Point Box<3>::physical( const Point &pos ) const
{
    return Point( d_range[0] + pos[0] * ( d_range[1] - d_range[0] ),
                  d_range[2] + pos[1] * ( d_range[3] - d_range[2] ),
                  d_range[4] + pos[2] * ( d_range[5] - d_range[4] ) );
}
template<std::size_t NDIM>
Point Grid<NDIM>::physical( const Point &pos ) const
{
    Point point( NDIM );
    for ( size_t d = 0; d < NDIM; d++ ) {
        double p = pos[d] * ( d_coord[d].size() - 1 );
        int i    = std::min<int>( p, d_coord[d].size() - 2 );
        point[d] = d_coord[d][i] + ( p - i ) * ( d_coord[d][i + 1] - d_coord[d][i] );
    }
    return point;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
template<>
Point Box<1>::logical( const Point &pos ) const
{
    return Point( ( pos[0] - d_range[0] ) / ( d_range[1] - d_range[0] ) );
}
template<>
Point Box<2>::logical( const Point &pos ) const
{
    return Point( ( pos[0] - d_range[0] ) / ( d_range[1] - d_range[0] ),
                  ( pos[1] - d_range[2] ) / ( d_range[3] - d_range[2] ) );
}
template<>
Point Box<3>::logical( const Point &pos ) const
{
    return Point( ( pos[0] - d_range[0] ) / ( d_range[1] - d_range[0] ),
                  ( pos[1] - d_range[2] ) / ( d_range[3] - d_range[2] ),
                  ( pos[2] - d_range[4] ) / ( d_range[5] - d_range[4] ) );
}
template<std::size_t NDIM>
Point Grid<NDIM>::logical( const Point &pos ) const
{
    Point logical;
    for ( size_t d = 0; d < NDIM; d++ ) {
        int i = AMP::Utilities::findfirst( d_coord[d], pos[d] );
        i     = std::max<int>( i, 1 );
        i     = std::min<int>( i, d_coord[d].size() - 1 );
        logical[d] =
            ( i - 1 ) + ( pos[d] - d_coord[d][i - 1] ) / ( d_coord[d][i] - d_coord[d][i - 1] );
        logical[d] /= d_coord[d].size();
    }
    return logical;
}


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
template<std::size_t NDIM>
void Box<NDIM>::displaceMesh( const double *x )
{
    for ( size_t d = 0; d < NDIM; d++ ) {
        d_range[2 * d + 0] += x[d];
        d_range[2 * d + 1] += x[d];
    }
}
template<std::size_t NDIM>
void Grid<NDIM>::displaceMesh( const double *x )
{
    for ( size_t d = 0; d < NDIM; d++ ) {
        d_range[2 * d + 0] += x[d];
        d_range[2 * d + 1] += x[d];
        for ( auto &tmp : d_coord[d] )
            tmp += x[d];
    }
}


/********************************************************
 * Explicit instantiations                               *
 ********************************************************/
template class Box<1>;
template class Box<2>;
template class Box<3>;
template class Grid<1>;
template class Grid<2>;
template class Grid<3>;


} // namespace Geometry
} // namespace AMP
