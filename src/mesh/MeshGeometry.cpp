#include "AMP/mesh/MeshGeometry.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshElement.h"
#include "AMP/mesh/MeshUtilities.h"
#include "AMP/mesh/shapes/GeometryHelpers.h"
#include "AMP/utils/kdtree2.h"


namespace AMP::Geometry {


/********************************************************
 * Constructors                                          *
 ********************************************************/
MeshGeometry::MeshGeometry( std::shared_ptr<AMP::Mesh::Mesh> mesh ) : d_mesh( std::move( mesh ) )
{
    AMP_ASSERT( d_mesh );
    AMP_ASSERT( static_cast<int>( d_mesh->getGeomType() ) == d_mesh->getDim() - 1 );
    d_physicalDim = d_mesh->getDim();
    // Get the block ids (they will translate to surface ids)
    d_surfaceIds = d_mesh->getBlockIDs();
    // Initialize position related data
    d_find = AMP::Mesh::ElementFinder( d_mesh );
}
std::unique_ptr<AMP::Geometry::Geometry> MeshGeometry::clone() const
{
    return std::unique_ptr<AMP::Geometry::Geometry>( new MeshGeometry( d_mesh->clone() ) );
}


/********************************************************
 * Get the distance to the surface                       *
 ********************************************************/
Point MeshGeometry::nearest( const Point &pos ) const { return d_find.nearest( pos ).second; }
double MeshGeometry::distance( const Point &pos, const Point &dir ) const
{
    return d_find.distance( pos, dir );
}
bool MeshGeometry::inside( const Point &pos ) const
{
    // Get the nearest element
    auto [elem, p] = d_find.nearest( pos );
    if ( p == pos )
        return true;
    //
    NULL_USE( elem );
    AMP_ERROR( "Not finished" );
    /*// If the nearest element is a surface, add the neighbors
    if ( elems.size() == 1 ) {
        auto neighbors = elems[0].getNeighbors();
        for ( const auto &tmp : neighbors )
            elems.push_back( *tmp );
    }
    // Search each element to determine if it is behind the plane
    bool test = true;
    for ( const auto &elem : elems ) {
        // Get the element normal
        auto n = elem.norm();
        // Get a vector to the surface
        auto vec = elem.centroid() - pos;
        // Check if the vector to intersection is in the same direction as the normal
        double t = dot( vec, n );
        test     = test && t >= -1e-12;
    }

    return test;*/
    return false;
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
    auto elem = d_find.nearest( x ).first;
    AMP_ASSERT( !elem.isNull() );
    for ( auto id : d_surfaceIds ) {
        if ( elem.isInBlock( id ) )
            return id;
    }
    return 0;
}
Point MeshGeometry::surfaceNorm( const Point &x ) const
{
    auto elem = d_find.nearest( x ).first;
    AMP_ASSERT( !elem.isNull() );
    return elem.norm();
}


/********************************************************
 * Get the centroid/box                                  *
 ********************************************************/
Point MeshGeometry::centroid() const
{
    // We are using the full mesh
    // Take the center of the box for now (this is not accurate)
    auto box = d_mesh->getBoundingBox();
    Point centroid( box.size() / 2 );
    for ( size_t d = 0; d < centroid.size(); d++ )
        centroid[d] = 0.5 * ( box[2 * d + 1] + box[2 * d] );
    return centroid;
}
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
 * Return the volume                                     *
 ********************************************************/
double MeshGeometry::volume() const
{
    std::vector<int> N( d_mesh->getDim(), 100 );
    auto V = AMP::Mesh::volumeOverlap( *this, N );
    return V.sum();
}


/********************************************************
 * Get the centroid/box                                  *
 ********************************************************/
void MeshGeometry::displace( const double *x )
{
    std::vector<double> x2( x, x + d_mesh->getDim() );
    d_mesh->displaceMesh( x2 );
}


/********************************************************
 * Get the nearest element                               *
 ********************************************************/
bool MeshGeometry::isConvex() const
{
    bool is_convex = true;
    auto [lb, ub]  = box();
    double tol     = 1e-4 * abs( ub - lb );
    for ( const auto &elem : d_mesh->getIterator( d_mesh->getGeomType() ) ) {
        // Get the normal to the plane of the element
        auto n = elem.norm();
        // Get a point in the plane
        auto a = elem.centroid();
        // Check the verticies of the neighboring elements to ensure they are not behind the plane
        for ( const auto &neighbor : elem.getNeighbors() ) {
            for ( const auto &node : neighbor->getElements( AMP::Mesh::GeomType::Vertex ) ) {
                auto p    = node.coord();
                double v  = dot( n, a - p );
                is_convex = is_convex && v >= -tol;
            }
        }
    }
    return is_convex;
}


/********************************************************
 * Compare the geometry                                  *
 ********************************************************/
bool MeshGeometry::operator==( const Geometry &rhs ) const
{
    auto geom = dynamic_cast<const MeshGeometry *>( &rhs );
    if ( !geom )
        return false;
    AMP_ERROR( "Not finished" );
    return false;
}


} // namespace AMP::Geometry
