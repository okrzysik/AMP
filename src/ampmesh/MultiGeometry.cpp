#include "AMP/ampmesh/MultiGeometry.h"
#include "AMP/ampmesh/Geometry.h"
#include "AMP/utils/UtilityMacros.h"


namespace AMP {
namespace Geometry {


MultiGeometry::MultiGeometry( const std::vector<Geometry::shared_ptr> &geom )
    : Geometry(), d_geom( geom )
{
    d_physicalDim = 0;
    for ( const auto &geom : d_geom )
        d_physicalDim = std::max( d_physicalDim, geom->getDim() );
}
Point MultiGeometry::nearest( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return {};
}
double MultiGeometry::distance( const Point &pos, const Point &dir ) const
{
    double dist = std::numeric_limits<double>::infinity();
    for ( const auto &geom : d_geom ) {
        double d  = geom->distance( pos, dir );
        double d1 = fabs( d );
        double d2 = fabs( dist );
        if ( fabs( d1 - d2 ) < 1e-6 * d1 ) {
            // The two distances are ~ equal, take the inside distance if it exists
            dist = std::min( d, dist );
        } else if ( d1 < d2 ) {
            // Keep the smallest distance
            dist = d;
        }
    }
    return dist;
}
bool MultiGeometry::inside( const Point &pos ) const
{
    bool inside = false;
    for ( const auto &geom : d_geom )
        inside = inside || geom->inside( pos );
    return inside;
}
int MultiGeometry::NSurface() const
{
    AMP_ERROR( "NSurface is not valid for MultiGeometry" );
    return 0;
}
int MultiGeometry::surface( const Point & ) const
{
    AMP_ERROR( "surface is not valid for MultiGeometry" );
    return 0;
}
Point MultiGeometry::surfaceNorm( const Point & ) const
{
    AMP_ERROR( "surfaceNorm is not valid for MultiGeometry" );
    return Point();
}
Point MultiGeometry::centroid() const
{
    auto range = box();
    return 0.5 * ( range.first + range.second );
}
std::pair<Point, Point> MultiGeometry::box() const
{
    Point lb = { 1e200, 1e200, 1e200 };
    Point ub = { -1e200, -1e200, -1e200 };
    for ( const auto &geom : d_geom ) {
        auto [lb2, ub2] = geom->box();
        lb.x()          = std::min( lb.x(), lb2.x() );
        lb.y()          = std::min( lb.y(), lb2.y() );
        lb.z()          = std::min( lb.z(), lb2.z() );
        ub.x()          = std::max( ub.x(), ub2.x() );
        ub.y()          = std::max( ub.y(), ub2.y() );
        ub.z()          = std::max( ub.z(), ub2.z() );
    }
    for ( int d = d_physicalDim; d < 3; d++ ) {
        lb[d] = 0;
        ub[d] = 0;
    }
    lb.setNdim( d_physicalDim );
    ub.setNdim( d_physicalDim );
    return std::make_pair( lb, ub );
}
double MultiGeometry::volume() const
{
    double V = 0;
    for ( const auto &geom : d_geom )
        V += geom->volume();
    return V;
}
void MultiGeometry::displace( const double *x )
{
    for ( const auto &geom : d_geom )
        geom->displace( x );
}
std::unique_ptr<AMP::Geometry::Geometry> MultiGeometry::clone() const
{
    std::vector<Geometry::shared_ptr> geom2;
    for ( const auto &geom : d_geom )
        geom2.push_back( geom->clone() );
    return std::make_unique<MultiGeometry>( geom2 );
}


/********************************************************
 * Compare the geometry                                  *
 ********************************************************/
bool MultiGeometry::operator==( const Geometry &rhs ) const
{
    auto geom = dynamic_cast<const MultiGeometry *>( &rhs );
    if ( !geom )
        return false;
    return d_geom == geom->d_geom;
}


} // namespace Geometry
} // namespace AMP
