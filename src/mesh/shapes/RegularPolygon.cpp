#include "AMP/mesh/shapes/RegularPolygon.h"
#include "AMP/mesh/shapes/GeometryHelpers.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
RegularPolygon::RegularPolygon( std::shared_ptr<const AMP::Database> db )
{
    d_physicalDim = 2;
    d_logicalDim  = 2;
    d_offset[0]   = 0;
    d_offset[1]   = 0;
    d_N           = db->getScalar<double>( "N" );
    d_R           = db->getScalar<double>( "R" );
    computeNorms();
}
RegularPolygon::RegularPolygon( int N, double R ) : LogicalGeometry(), d_N( N ), d_R( R )
{
    d_physicalDim = 2;
    d_logicalDim  = 2;
    d_offset[0]   = 0;
    d_offset[1]   = 0;
    computeNorms();
}
void RegularPolygon::computeNorms()
{
    // Get the verticies
    d_verticies = GeometryHelpers::get_poly_verticies( d_N, d_R );
    for ( auto &p : d_verticies ) {
        p[0] += d_offset[0];
        p[1] += d_offset[1];
    }
    // Calculate the normals
    auto p0       = centroid();
    auto calcNorm = [p0]( const Point &v1, const Point &v2 ) {
        auto v  = normalize( v2 - v1 );
        Point n = { v.x(), -v.y() };
        if ( dot( n, v2 - p0 ) < 0 )
            return -n;
        return n;
    };
    d_norm.resize( d_N );
    d_norm[0] = calcNorm( d_verticies.back(), d_verticies[0] );
    for ( size_t i = 1; i < d_verticies.size(); i++ )
        d_norm[i] = calcNorm( d_verticies[i - 1], d_verticies[i] );
}


/********************************************************
 * Compute the nearest point on the surface              *
 ********************************************************/
std::tuple<Point, double, int> RegularPolygon::nearest2( const Point &pos ) const
{
    std::array<double, 2> p0 = { pos.x(), pos.y() };
    // Check the intersection with each line segment (keeping the closest)
    Point p  = GeometryHelpers::nearest( d_verticies[0], d_verticies.back(), p0 );
    double d = ( p - pos ).norm();
    int k    = 0;
    for ( size_t i = 1; i < d_verticies.size(); i++ ) {
        Point p2  = GeometryHelpers::nearest( d_verticies[i], d_verticies[i - 1], p0 );
        double d2 = ( p2 - pos ).norm();
        if ( d2 < d ) {
            d = d2;
            p = p2;
            k = i;
        }
    }
    return std::tie( p, d, k );
}
Point RegularPolygon::nearest( const Point &pos ) const
{
    if ( inside( pos ) )
        return pos;
    auto r = nearest2( pos );
    return std::get<0>( r );
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double RegularPolygon::distance( const Point &pos, const Point &ang ) const
{
    // Check the intersection with each line segment (keeping the closest)
    double d = GeometryHelpers::distanceToLine( pos, ang, d_verticies[0], d_verticies.back() );
    for ( size_t i = 1; i < d_verticies.size(); i++ ) {
        double d2 = GeometryHelpers::distanceToLine( pos, ang, d_verticies[i], d_verticies[i - 1] );
        if ( d2 < d )
            d = d2;
    }
    if ( d == std::numeric_limits<double>::infinity() )
        return d;
    if ( inside( pos ) )
        return -d;
    return d;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool RegularPolygon::inside( const Point &pos ) const
{
    auto L    = logical( pos );
    double t1 = -1e-12;
    double t2 = 1.0 + 1e-12;
    return L.x() >= t1 && L.y() >= t1 && L.z() >= t1 && L.x() <= t2 && L.y() <= t2 && L.z() <= t2;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
int RegularPolygon::NSurface() const { return d_N; }
int RegularPolygon::surface( const Point &pos ) const
{
    auto r = nearest2( pos );
    return std::get<2>( r );
}
Point RegularPolygon::surfaceNorm( const Point &pos ) const
{
    int i = surface( pos );
    return d_norm[i];
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point RegularPolygon::physical( const Point &pos ) const
{
    auto tmp = GeometryHelpers::map_logical_poly( d_N, d_R, pos.x(), pos.y() );
    double x = tmp.first + d_offset[0];
    double y = tmp.second + d_offset[1];
    return { x, y };
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point RegularPolygon::logical( const Point &pos ) const
{
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    auto tmp = GeometryHelpers::map_poly_logical( d_N, d_R, x, y );
    return Point( tmp.first, tmp.second );
}


/********************************************************
 * Return the centroid and bounding box                  *
 ********************************************************/
Point RegularPolygon::centroid() const { return { d_offset[0], d_offset[1] }; }
std::pair<Point, Point> RegularPolygon::box() const
{
    Point lb = { d_offset[0] - d_R, d_offset[1] - d_R };
    Point ub = { d_offset[0] + d_R, d_offset[1] + d_R };
    return { lb, ub };
}


/********************************************************
 * Return the volume                                     *
 ********************************************************/
double RegularPolygon::volume() const
{
    constexpr double pi = 3.141592653589793;
    return 0.5 * d_N * d_R * d_R * sin( 2 * pi / d_N );
}


/********************************************************
 * Return the logical grid                               *
 ********************************************************/
std::vector<int> RegularPolygon::getLogicalGridSize( const std::vector<int> &x ) const
{
    AMP_INSIST( x.size() == 1u, "Size must be an array of length 1" );
    return { 2 * x[0], 2 * x[0] };
}
std::vector<int> RegularPolygon::getLogicalGridSize( const std::vector<double> &res ) const
{
    AMP_INSIST( res.size() == 2u, "Resolution must be an array of length 2" );
    return { (int) ( d_R / res[0] ), (int) ( d_R / res[1] ) };
}
std::vector<bool> RegularPolygon::getPeriodicDim() const { return { false, false }; }
std::vector<int> RegularPolygon::getLogicalSurfaceIds() const { return { 1, 1, 1, 1 }; }


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void RegularPolygon::displace( const double *x )
{
    // Update the offsets
    d_offset[0] += x[0];
    d_offset[1] += x[1];
    // Get the offsets
    d_verticies = GeometryHelpers::get_poly_verticies( d_N, d_R );
    for ( auto &p : d_verticies ) {
        p[0] += d_offset[0];
        p[1] += d_offset[1];
    }
}


/********************************************************
 * Clone the object                                      *
 ********************************************************/
std::unique_ptr<AMP::Geometry::Geometry> RegularPolygon::clone() const
{
    return std::make_unique<RegularPolygon>( *this );
}


/********************************************************
 * Compare the geometry                                  *
 ********************************************************/
bool RegularPolygon::operator==( const Geometry &rhs ) const
{
    auto geom = dynamic_cast<const RegularPolygon *>( &rhs );
    if ( !geom )
        return false;
    return d_N == geom->d_N && d_R == geom->d_R && d_offset == geom->d_offset &&
           d_verticies == geom->d_verticies && d_norm == geom->d_norm;
}


} // namespace Geometry
} // namespace AMP
