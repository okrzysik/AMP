#include "AMP/ampmesh/MultiGeometry.h"


namespace AMP {
namespace Geometry {


MultiGeometry::MultiGeometry( const std::vector<Geometry::shared_ptr> &geom ) : d_geom( geom ) {}


uint8_t MultiGeometry::getDim() const
{
    uint8_t dim = 0;
    for ( const auto &geom : d_geom )
        dim = std::max( dim, geom->getDim() );
    return dim;
}
double MultiGeometry::distance( const Point &pos, const Point &dir ) const
{
    double dist = std::numeric_limits<double>::infinity();
    for ( const auto &geom : d_geom )
        dist = std::min( dist, geom->distance( pos, dir ) );
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
void MultiGeometry::displaceMesh( const double *x )
{
    for ( const auto &geom : d_geom )
        geom->displaceMesh( x );
}


} // namespace Geometry
} // namespace AMP
