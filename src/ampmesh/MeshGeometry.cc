#include "AMP/ampmesh/MeshGeometry.h"
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshElement.h"


namespace AMP::Geometry {


/********************************************************
 * Constructors                                          *
 ********************************************************/
MeshGeometry::MeshGeometry( std::unique_ptr<AMP::Mesh::Mesh> mesh ) : d_mesh( std::move( mesh ) )
{
    AMP_ASSERT( d_mesh );
    AMP_ASSERT( static_cast<int>( d_mesh->getGeomType() ) == getDim() - 1 );
    d_physicalDim = d_mesh->getDim();
    initialize();
}
std::unique_ptr<AMP::Geometry::Geometry> MeshGeometry::clone() const
{
    return std::unique_ptr<AMP::Geometry::Geometry>( new MeshGeometry( d_mesh->clone() ) );
}
void MeshGeometry::initialize()
{
    // Compute the centroid
    // Take the center of the box for now (this is not accurate)
    auto b     = box();
    d_centroid = Point( b.first.size() );
    for ( size_t d = 0; d < d_centroid.size(); d++ )
        d_centroid[d] = 0.5 * ( b.first[d] + b.second[d] );
    // Get the block ids (they will translate to surface ids)
    d_surfaceIds = d_mesh->getBlockIDs();
}


/********************************************************
 * Get the distance to the surface                       *
 ********************************************************/
Point MeshGeometry::nearest( const Point &pos ) const
{
    NULL_USE( pos );
    AMP_ERROR( "Not finished" );
    return {};
}
double MeshGeometry::distance( const Point &pos, const Point &dir ) const
{
    NULL_USE( pos );
    NULL_USE( dir );
    AMP_ERROR( "Not finished" );
    return 0;
}
bool MeshGeometry::inside( const Point &pos ) const
{
    double dist = distance( pos, pos - d_centroid );
    return dist <= 0;
}


/********************************************************
 * Get the surface                                       *
 ********************************************************/
int MeshGeometry::NSurface() const { return d_surfaceIds.size(); }
int MeshGeometry::surface( const Point &x ) const
{
    if ( d_surfaceIds.empty() )
        return 0;
    if ( d_surfaceIds.size() == 1 )
        return d_surfaceIds[0];
    auto elem = getNearest( x );
    for ( auto id : d_surfaceIds ) {
        if ( elem.isInBlock( id ) )
            return id;
    }
    return 0;
}
Point MeshGeometry::surfaceNorm( const Point &x ) const
{
    auto elem = getNearest( x );
    return elem.norm();
}


/********************************************************
 * Get the centroid/box                                  *
 ********************************************************/
Point MeshGeometry::centroid() const { return d_centroid; }
std::pair<Point, Point> MeshGeometry::box() const
{
    auto box = d_mesh->getBoundingBox();
    Point p0( box.size() / 2 ), p1( box.size() / 2 );
    for ( size_t d = 0; d < box.size() / 2; d++ ) {
        p0[d] = box[2 * d];
        p1[d] = box[2 * d + 1];
    }
    return std::make_pair( p0, p1 );
}


/********************************************************
 * Get the centroid/box                                  *
 ********************************************************/
void MeshGeometry::displace( const double *x )
{
    std::vector<double> x2( d_mesh->getDim() );
    for ( size_t i = 0; i < x2.size(); i++ )
        x2[i] = x[i];
    d_mesh->displaceMesh( x2 );
}


/********************************************************
 * Get the nearest element                               *
 ********************************************************/
AMP::Mesh::MeshElement MeshGeometry::getNearest( const Point &x ) const
{
    NULL_USE( x );
    AMP_ERROR( "Not finished" );
    return AMP::Mesh::MeshElement();
}


/********************************************************
 * Get the nearest element                               *
 ********************************************************/
bool MeshGeometry::isConvex() const
{
    AMP_ERROR( "Not finished" );
    return false;
}


} // namespace AMP::Geometry
