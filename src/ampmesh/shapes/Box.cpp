#include "AMP/ampmesh/shapes/Box.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"


#include <algorithm>
#include <array>
#include <vector>


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructors                                          *
 ********************************************************/
template<std::size_t NDIM>
Box<NDIM>::Box() : LogicalGeometry()
{
    static_assert( NDIM >= 0 && NDIM <= 3, "Invalid number of dimensions" );
    d_physicalDim = NDIM;
    d_logicalDim  = NDIM;
    d_range.fill( 0 );
}
template<std::size_t NDIM>
Box<NDIM>::Box( std::shared_ptr<const AMP::Database> db ) : Box<NDIM>()
{
    auto range = db->getVector<double>( "Range" );
    AMP_INSIST( range.size() == 2 * NDIM, "Range must be 2*dim for cube generator" );
    for ( size_t i = 0; i < 2 * NDIM; i++ )
        d_range[i] = range[i];
}
template<std::size_t NDIM>
Box<NDIM>::Box( const std::vector<double> &range ) : Box<NDIM>()
{
    AMP_INSIST( range.size() == 2 * NDIM, "Range must be 2*dim" );
    for ( size_t i = 0; i < 2 * NDIM; i++ )
        d_range[i] = range[i];
}
template<std::size_t NDIM>
Grid<NDIM>::Grid( std::shared_ptr<const AMP::Database> db ) : Box<NDIM>()
{
    const char *names[3] = { "x_grid", "y_grid", "z_grid" };
    for ( size_t d = 0; d < NDIM; d++ ) {
        this->d_coord[d] = std::move( db->getVector<double>( names[d] ) );
        AMP_ASSERT( d_coord[d].size() >= 2 );
        for ( size_t j = 1; j < d_coord[d].size(); j++ )
            AMP_ASSERT( d_coord[d][j] > d_coord[d][j - 1] );
        this->d_range[2 * d + 0] = 1e100;
        this->d_range[2 * d + 1] = -1e100;
        for ( const auto tmp : d_coord[d] ) {
            this->d_range[2 * d + 0] = std::min( this->d_range[2 * d + 0], tmp );
            this->d_range[2 * d + 1] = std::max( this->d_range[2 * d + 1], tmp );
        }
    }
}
template<std::size_t NDIM>
Grid<NDIM>::Grid( const std::vector<std::vector<double>> &coord ) : Box<NDIM>()
{
    AMP_ASSERT( coord.size() == NDIM );
    for ( size_t d = 0; d < NDIM; d++ ) {
        d_coord[d] = coord[d];
        AMP_ASSERT( d_coord[d].size() > 2 );
        for ( size_t j = 1; j < d_coord[d].size(); j++ )
            AMP_ASSERT( d_coord[d][j] > d_coord[d][j - 1] );
        this->d_range[2 * d + 0] = 1e100;
        this->d_range[2 * d + 1] = -1e100;
        for ( const auto tmp : d_coord[d] ) {
            this->d_range[2 * d + 0] = std::min( this->d_range[2 * d + 0], tmp );
            this->d_range[2 * d + 1] = std::max( this->d_range[2 * d + 1], tmp );
        }
    }
}


/********************************************************
 * Get the object name                                   *
 ********************************************************/
template<std::size_t NDIM>
std::string Box<NDIM>::getName() const
{
    if constexpr ( NDIM == 1 )
        return "Box<1>";
    else if constexpr ( NDIM == 2 )
        return "Box<2>";
    else if constexpr ( NDIM == 3 )
        return "Box<3>";
}
template<std::size_t NDIM>
std::string Grid<NDIM>::getName() const
{
    if constexpr ( NDIM == 1 )
        return "Grid<1>";
    else if constexpr ( NDIM == 2 )
        return "Grid<2>";
    else if constexpr ( NDIM == 3 )
        return "Grid<3>";
}


/********************************************************
 * Compute the nearest point on the surface              *
 ********************************************************/
template<std::size_t NDIM>
Point Box<NDIM>::nearest( const Point &pos ) const
{
    auto p = logical( pos );
    p.x()  = std::min( p.x(), 1.0 );
    p.x()  = std::max( p.x(), 0.0 );
    p.y()  = std::min( p.y(), 1.0 );
    p.y()  = std::max( p.y(), 0.0 );
    p.z()  = std::min( p.z(), 1.0 );
    p.z()  = std::max( p.z(), 0.0 );
    return physical( p );
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
template<std::size_t NDIM>
double Box<NDIM>::distance( const Point &pos, const Point &ang ) const
{
    return GeometryHelpers::distanceToBox( pos, ang, d_range );
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
    for ( size_t i = 1; i < 2 * NDIM; i++ ) {
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
    Point logical( NDIM );
    for ( size_t d = 0; d < NDIM; d++ ) {
        int i = AMP::Utilities::findfirst( d_coord[d], pos[d] );
        i     = std::max<int>( i, 1 );
        i     = std::min<int>( i, d_coord[d].size() - 1 );
        logical[d] =
            ( i - 1 ) + ( pos[d] - d_coord[d][i - 1] ) / ( d_coord[d][i] - d_coord[d][i - 1] );
        logical[d] /= ( d_coord[d].size() - 1 );
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
std::pair<Point, Point> Box<NDIM>::box() const
{
    Point lb( NDIM, { d_range[0], d_range[2], d_range[4] } );
    Point ub( NDIM, { d_range[1], d_range[3], d_range[5] } );
    return { lb, ub };
}


/********************************************************
 * Return the volume                                     *
 ********************************************************/
template<std::size_t NDIM>
double Box<NDIM>::volume() const
{
    double v = 1.0;
    for ( size_t d = 0; d < NDIM; d++ )
        v *= d_range[2 * d + 1] - d_range[2 * d];
    return v;
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
std::vector<int> Box<NDIM>::getLogicalGridSize( const std::vector<double> &res ) const
{
    AMP_ASSERT( res.size() == NDIM );
    std::vector<int> size( NDIM );
    for ( size_t d = 0; d < NDIM; d++ )
        size[d] = ( d_range[2 * d + 1] - d_range[2 * d] ) / res[d];
    return size;
}
template<std::size_t NDIM>
std::vector<bool> Box<NDIM>::getPeriodicDim() const
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


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
template<std::size_t NDIM>
void Box<NDIM>::displace( const double *x )
{
    for ( size_t d = 0; d < NDIM; d++ ) {
        d_range[2 * d + 0] += x[d];
        d_range[2 * d + 1] += x[d];
    }
}
template<std::size_t NDIM>
void Grid<NDIM>::displace( const double *x )
{
    Box<NDIM>::displace( x );
    for ( size_t d = 0; d < NDIM; d++ ) {
        for ( auto &tmp : d_coord[d] )
            tmp += x[d];
    }
}


/********************************************************
 * Clone the object                                      *
 ********************************************************/
template<std::size_t NDIM>
std::unique_ptr<AMP::Geometry::Geometry> Box<NDIM>::clone() const
{
    return std::make_unique<Box<NDIM>>( *this );
}
template<std::size_t NDIM>
std::unique_ptr<AMP::Geometry::Geometry> Grid<NDIM>::clone() const
{
    return std::make_unique<Grid<NDIM>>( *this );
}


/********************************************************
 * Compare the geometry                                  *
 ********************************************************/
template<std::size_t NDIM>
bool Box<NDIM>::operator==( const Geometry &rhs ) const
{
    auto geom = dynamic_cast<const Box<NDIM> *>( &rhs );
    if ( !geom )
        return false;
    return d_range == geom->d_range;
}
template<std::size_t NDIM>
bool Grid<NDIM>::operator==( const Geometry &rhs ) const
{
    auto geom = dynamic_cast<const Grid<NDIM> *>( &rhs );
    if ( !geom )
        return false;
    bool test = this->d_range == geom->d_range;
    for ( size_t d = 0; d < NDIM; d++ )
        test = test && d_coord[d] == geom->d_coord[d];
    return test;
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
