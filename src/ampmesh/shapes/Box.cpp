#include "AMP/ampmesh/shapes/Box.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"


#include <algorithm>
#include <vector>


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructors                                          *
 ********************************************************/
template<std::size_t NDIM>
Box<NDIM>::Box( AMP::shared_ptr<AMP::Database> db )
{
    d_physicalDim = NDIM;
    d_logicalDim  = NDIM;
    static_assert( NDIM >= 0 && NDIM <= 3, "Invalid number of dimensions" );
    auto range = db->getDoubleArray( "Range" );
    AMP_INSIST( range.size() == 2 * NDIM, "Range must be 2*dim for cube generator" );
    for ( int i = 0; i < 3; i++ ) {
        d_range[2 * i + 0] = -1e100;
        d_range[2 * i + 1] = 1e100;
    }
    for ( size_t i = 0; i < 2 * NDIM; i++ )
        d_range[i] = range[i];
}
template<std::size_t NDIM>
Box<NDIM>::Box( const std::vector<double> &range ) : Geometry()
{
    d_physicalDim = NDIM;
    d_logicalDim  = NDIM;
    static_assert( NDIM >= 0 && NDIM <= 3, "Invalid number of dimensions" );
    for ( int i = 0; i < 3; i++ ) {
        d_range[2 * i + 0] = -1e100;
        d_range[2 * i + 1] = 1e100;
    }
    for ( size_t i = 0; i < 2 * NDIM; i++ )
        d_range[i] = range[i];
}
template<std::size_t NDIM>
Grid<NDIM>::Grid( AMP::shared_ptr<AMP::Database> db ) : Geometry()
{
    d_physicalDim = NDIM;
    d_logicalDim  = NDIM;
    static_assert( NDIM >= 0 && NDIM <= 3, "Invalid number of dimensions" );
    for ( int i = 0; i < 3; i++ ) {
        d_range[2 * i + 0] = -1e100;
        d_range[2 * i + 1] = 1e100;
    }
    const char *names[3] = { "x_grid", "y_grid", "z_grid" };
    for ( size_t i = 0; i < NDIM; i++ ) {
        d_coord[i] = std::move( db->getDoubleArray( names[i] ) );
        for ( const auto tmp : d_coord[i] ) {
            d_range[2 * i + 0] = std::min( d_range[2 * i + 0], tmp );
            d_range[2 * i + 1] = std::max( d_range[2 * i + 1], tmp );
        }
    }
}
template<std::size_t NDIM>
Grid<NDIM>::Grid( const std::vector<std::vector<double>> &coord ) : Geometry()
{
    d_physicalDim = NDIM;
    d_logicalDim  = NDIM;
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
 * Get the object name                                   *
 ********************************************************/
template<>
std::string Box<1>::getName() const
{
    return "Box<1>";
}
template<>
std::string Box<2>::getName() const
{
    return "Box<2>";
}
template<>
std::string Box<3>::getName() const
{
    return "Box<3>";
}
template<>
std::string Grid<1>::getName() const
{
    return "Grid<1>";
}
template<>
std::string Grid<2>::getName() const
{
    return "Grid<2>";
}
template<>
std::string Grid<3>::getName() const
{
    return "Grid<3>";
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
static inline double calcDist( const Point &pos, const Point &ang, const double *range )
{
    constexpr double tol = 1e-12;
    // Compute the distance to each surface
    double d1 = ( range[0] - pos.x() ) / ang.x();
    double d2 = ( range[1] - pos.x() ) / ang.x();
    double d3 = ( range[2] - pos.y() ) / ang.y();
    double d4 = ( range[3] - pos.y() ) / ang.y();
    double d5 = ( range[4] - pos.z() ) / ang.z();
    double d6 = ( range[5] - pos.z() ) / ang.z();
    if ( d1 < 0 )
        d1 = std::numeric_limits<double>::infinity();
    if ( d2 < 0 )
        d2 = std::numeric_limits<double>::infinity();
    if ( d3 < 0 )
        d3 = std::numeric_limits<double>::infinity();
    if ( d4 < 0 )
        d4 = std::numeric_limits<double>::infinity();
    if ( d5 < 0 )
        d5 = std::numeric_limits<double>::infinity();
    if ( d6 < 0 )
        d6 = std::numeric_limits<double>::infinity();
    // Check if the intersection of each surface is within the bounds of the box
    auto p1     = pos + d1 * ang;
    auto p2     = pos + d2 * ang;
    auto p3     = pos + d3 * ang;
    auto p4     = pos + d4 * ang;
    auto p5     = pos + d5 * ang;
    auto p6     = pos + d6 * ang;
    auto inside = [tol, range]( const Point &p ) {
        return ( ( p.x() >= range[0] - tol ) && ( p.x() <= range[1] + tol ) ) &&
               ( ( p.y() >= range[2] - tol ) && ( p.y() <= range[3] + tol ) ) &&
               ( ( p.z() >= range[4] - tol ) && ( p.z() <= range[5] + tol ) );
    };
    if ( !inside( p1 ) )
        d1 = std::numeric_limits<double>::infinity();
    if ( !inside( p2 ) )
        d2 = std::numeric_limits<double>::infinity();
    if ( !inside( p3 ) )
        d3 = std::numeric_limits<double>::infinity();
    if ( !inside( p4 ) )
        d4 = std::numeric_limits<double>::infinity();
    if ( !inside( p5 ) )
        d5 = std::numeric_limits<double>::infinity();
    if ( !inside( p6 ) )
        d6 = std::numeric_limits<double>::infinity();
    // Return the closest surface
    double d = std::min( { d1, d2, d3, d4, d5, d6 } );
    if ( inside( pos ) && d < 1e100 )
        d = -d;
    return d;
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
    double x  = ( pos.x() - d_range[0] ) / ( d_range[1] - d_range[0] );
    double t1 = -1e-12;
    double t2 = 1.0 + 1e-12;
    return x >= t1 && x <= t2;
}
template<>
bool Box<2>::inside( const Point &pos ) const
{
    double x  = ( pos.x() - d_range[0] ) / ( d_range[1] - d_range[0] );
    double y  = ( pos.y() - d_range[2] ) / ( d_range[3] - d_range[2] );
    double t1 = -1e-12;
    double t2 = 1.0 + 1e-12;
    return x >= t1 && x <= t2 && y >= t1 && y <= t2;
}
template<>
bool Box<3>::inside( const Point &pos ) const
{
    double x  = ( pos.x() - d_range[0] ) / ( d_range[1] - d_range[0] );
    double y  = ( pos.y() - d_range[2] ) / ( d_range[3] - d_range[2] );
    double z  = ( pos.z() - d_range[4] ) / ( d_range[5] - d_range[4] );
    double t1 = -1e-12;
    double t2 = 1.0 + 1e-12;
    return x >= t1 && x <= t2 && y >= t1 && y <= t2 && z >= t1 && z <= t2;
}
template<>
bool Grid<1>::inside( const Point &pos ) const
{
    double x  = ( pos.x() - d_range[0] ) / ( d_range[1] - d_range[0] );
    double t1 = -1e-12;
    double t2 = 1.0 + 1e-12;
    return x >= t1 && x <= t2;
}
template<>
bool Grid<2>::inside( const Point &pos ) const
{
    double x  = ( pos.x() - d_range[0] ) / ( d_range[1] - d_range[0] );
    double y  = ( pos.y() - d_range[2] ) / ( d_range[3] - d_range[2] );
    double t1 = -1e-12;
    double t2 = 1.0 + 1e-12;
    return x >= t1 && x <= t2 && y >= t1 && y <= t2;
}
template<>
bool Grid<3>::inside( const Point &pos ) const
{
    double x  = ( pos.x() - d_range[0] ) / ( d_range[1] - d_range[0] );
    double y  = ( pos.y() - d_range[2] ) / ( d_range[3] - d_range[2] );
    double z  = ( pos.z() - d_range[4] ) / ( d_range[5] - d_range[4] );
    double t1 = -1e-12;
    double t2 = 1.0 + 1e-12;
    return x >= t1 && x <= t2 && y >= t1 && y <= t2 && z >= t1 && z <= t2;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
template<std::size_t NDIM>
int Box<NDIM>::surface( const Point &pos ) const
{
    double d[6] = { fabs( pos.x() - d_range[0] ), fabs( pos.x() - d_range[1] ),
                    fabs( pos.y() - d_range[2] ), fabs( pos.y() - d_range[3] ),
                    fabs( pos.z() - d_range[4] ), fabs( pos.z() - d_range[5] ) };
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
    double d[6] = { fabs( pos.x() - d_range[0] ), fabs( pos.x() - d_range[1] ),
                    fabs( pos.y() - d_range[2] ), fabs( pos.y() - d_range[3] ),
                    fabs( pos.z() - d_range[4] ), fabs( pos.z() - d_range[5] ) };
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
    return Point( d_range[0] + pos.x() * ( d_range[1] - d_range[0] ) );
}
template<>
Point Box<2>::physical( const Point &pos ) const
{
    return Point( d_range[0] + pos.x() * ( d_range[1] - d_range[0] ),
                  d_range[2] + pos.y() * ( d_range[3] - d_range[2] ) );
}
template<>
Point Box<3>::physical( const Point &pos ) const
{
    return Point( d_range[0] + pos.x() * ( d_range[1] - d_range[0] ),
                  d_range[2] + pos.y() * ( d_range[3] - d_range[2] ),
                  d_range[4] + pos.z() * ( d_range[5] - d_range[4] ) );
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
    return Point( ( pos.x() - d_range[0] ) / ( d_range[1] - d_range[0] ) );
}
template<>
Point Box<2>::logical( const Point &pos ) const
{
    return Point( ( pos.x() - d_range[0] ) / ( d_range[1] - d_range[0] ),
                  ( pos.y() - d_range[2] ) / ( d_range[3] - d_range[2] ) );
}
template<>
Point Box<3>::logical( const Point &pos ) const
{
    return Point( ( pos.x() - d_range[0] ) / ( d_range[1] - d_range[0] ),
                  ( pos.y() - d_range[2] ) / ( d_range[3] - d_range[2] ),
                  ( pos.z() - d_range[4] ) / ( d_range[5] - d_range[4] ) );
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
 * Return the centroid and bounding box                  *
 ********************************************************/
template<std::size_t NDIM>
Point Box<NDIM>::centroid() const
{
    Point p( NDIM );
    p.x() = 0.5 * ( d_range[0] + d_range[1] );
    p.y() = 0.5 * ( d_range[2] + d_range[3] );
    p.z() = 0.5 * ( d_range[4] + d_range[5] );
    return p;
}
template<std::size_t NDIM>
Point Grid<NDIM>::centroid() const
{
    Point p( NDIM );
    p.x() = 0.5 * ( d_range[0] + d_range[1] );
    p.y() = 0.5 * ( d_range[2] + d_range[3] );
    p.z() = 0.5 * ( d_range[4] + d_range[5] );
    return p;
}
template<std::size_t NDIM>
std::pair<Point, Point> Box<NDIM>::box() const
{
    Point lb( NDIM, { d_range[0], d_range[2], d_range[4] } );
    Point ub( NDIM, { d_range[1], d_range[3], d_range[5] } );
    return { lb, ub };
}
template<std::size_t NDIM>
std::pair<Point, Point> Grid<NDIM>::box() const
{
    Point lb( NDIM, { d_range[0], d_range[2], d_range[4] } );
    Point ub( NDIM, { d_range[1], d_range[3], d_range[5] } );
    return { lb, ub };
}


/********************************************************
 * Return the logical grid                               *
 ********************************************************/
template<std::size_t NDIM>
std::vector<int> Box<NDIM>::getLogicalGridSize( const std::vector<int> &x ) const
{
    AMP_ASSERT( x.size() == NDIM );
    return x;
}
template<std::size_t NDIM>
std::vector<int> Grid<NDIM>::getLogicalGridSize( const std::vector<int> &x ) const
{
    AMP_ASSERT( x.size() == NDIM );
    return x;
}
template<std::size_t NDIM>
std::vector<bool> Box<NDIM>::getPeriodicDim() const
{
    return std::vector<bool>( NDIM, false );
}
template<std::size_t NDIM>
std::vector<bool> Grid<NDIM>::getPeriodicDim() const
{
    return std::vector<bool>( NDIM, false );
}
template<std::size_t NDIM>
std::vector<int> Box<NDIM>::getLogicalSurfaceIds() const
{
    std::vector<int> ids( 2 * NDIM );
    for ( size_t i = 0; i < ids.size(); i++ )
        ids[i] = i;
    return ids;
}
template<std::size_t NDIM>
std::vector<int> Grid<NDIM>::getLogicalSurfaceIds() const
{
    std::vector<int> ids( 2 * NDIM );
    for ( size_t i = 0; i < ids.size(); i++ )
        ids[i] = i;
    return ids;
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
 * Clone the object                                      *
 ********************************************************/
template<std::size_t NDIM>
AMP::shared_ptr<AMP::Geometry::Geometry> Box<NDIM>::clone() const
{
    return AMP::make_shared<Box<NDIM>>( *this );
}
template<std::size_t NDIM>
AMP::shared_ptr<AMP::Geometry::Geometry> Grid<NDIM>::clone() const
{
    return AMP::make_shared<Grid<NDIM>>( *this );
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
