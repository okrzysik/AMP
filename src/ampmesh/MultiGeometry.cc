#include "AMP/ampmesh/MultiGeometry.h"
#include "AMP/ampmesh/Geometry.h"


namespace AMP {
namespace Geometry {


MultiGeometry::MultiGeometry( const std::vector<Geometry::shared_ptr> &geom )
    : Geometry(), d_geom( geom )
{
    d_physicalDim = 0;
    for ( const auto &geom : d_geom )
        d_physicalDim = std::max( d_physicalDim, geom->getDim() );
    d_logicalDim = 0;
}
double MultiGeometry::distance( const Point &pos, const Point &dir ) const
{
    double dist = std::numeric_limits<double>::infinity();
    for ( const auto &geom : d_geom ) {
        double d = geom->distance( pos, dir );
        if ( fabs( d ) < fabs( dist ) )
            dist = d;
        if ( d < 1e100 ) {
            bool test = geom->inside( pos + fabs( d ) * dir );
            geom->distance( pos, dir );
            geom->inside( pos + fabs( d ) * dir );
            AMP_ASSERT( test );
        }
    }
    if ( dist < 1e100 ) {
        bool test = inside( pos + fabs( dist ) * dir );
        inside( pos + fabs( dist ) * dir );
        AMP_ASSERT( test );
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
Point MultiGeometry::logical( const Point & ) const
{
    AMP_ERROR( "Converting to logical coordinates is not valid for MultiGeometry" );
    return Point();
}
Point MultiGeometry::physical( const Point & ) const
{
    AMP_ERROR( "Converting to physical coordinates is not valid for MultiGeometry" );
    return Point();
}
Point MultiGeometry::centroid() const
{
    auto range = box();
    return 0.5 * ( range.first + range.second );
}
std::pair<Point, Point> MultiGeometry::box() const
{
    std::pair<Point, Point> range;
    range.first  = { 1e200, 1e200, 1e200 };
    range.second = { -1e200, -1e200, -1e200 };
    for ( const auto &geom : d_geom ) {
        auto range2      = geom->box();
        range.first.x()  = std::min( range.first.x(), range2.first.x() );
        range.first.y()  = std::min( range.first.y(), range2.first.y() );
        range.first.z()  = std::min( range.first.z(), range2.first.z() );
        range.second.x() = std::max( range.second.x(), range2.second.x() );
        range.second.y() = std::max( range.second.y(), range2.second.y() );
        range.second.z() = std::max( range.second.z(), range2.second.z() );
    }
    for ( int d = d_physicalDim; d < 3; d++ ) {
        range.first[d]  = 0;
        range.second[d] = 0;
    }
    range.first.setNdim( d_physicalDim );
    range.second.setNdim( d_physicalDim );
    return range;
}
void MultiGeometry::displaceMesh( const double *x )
{
    for ( const auto &geom : d_geom )
        geom->displaceMesh( x );
}


} // namespace Geometry
} // namespace AMP
