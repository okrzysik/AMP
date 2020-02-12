#include "AMP/ampmesh/MeshGeometry.h"
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshElement.h"


namespace AMP::Geometry {


/********************************************************
 * Constructors                                          *
 ********************************************************/
MeshGeometry::MeshGeometry( std::shared_ptr<AMP::Mesh::Mesh> mesh ) : d_mesh( std::move( mesh ) )
{
    AMP_ASSERT( d_mesh );
    AMP_ASSERT( static_cast<int>( d_mesh->getGeomType() ) == d_mesh->getDim() - 1 );
    d_physicalDim = d_mesh->getDim();
    initialize();
}
std::unique_ptr<AMP::Geometry::Geometry> MeshGeometry::clone() const
{
    return std::unique_ptr<AMP::Geometry::Geometry>( new MeshGeometry( d_mesh->clone() ) );
}
void MeshGeometry::initialize()
{
    // Get the block ids (they will translate to surface ids)
    d_surfaceIds = d_mesh->getBlockIDs();
    // Get the nodes
    auto it = d_mesh->getIterator( AMP::Mesh::GeomType::Vertex );
    d_nodes.resize( 0 );
    d_nodes.reserve( it.size() );
    for ( auto &elem : it )
        d_nodes.push_back( elem );
    // Initialize position related data
    initializePosition();
}
void MeshGeometry::initializePosition() const
{
    // Update the internal hash
    d_pos_hash = d_mesh->positionHash();
    // Compute the centroid
    // Take the center of the box for now (this is not accurate)
    auto b     = box();
    d_centroid = Point( b.first.size() );
    for ( size_t d = 0; d < d_centroid.size(); d++ )
        d_centroid[d] = 0.5 * ( b.first[d] + b.second[d] );
    // Create a kdtree with the node positions
    std::vector<Point> pos( d_nodes.size() );
    for ( size_t i = 0; i < d_nodes.size(); i++ )
        pos[i] = d_nodes[i].coord();
    d_tree = kdtree( pos );
}

/********************************************************
 * Get the distance to the surface                       *
 ********************************************************/
Point MeshGeometry::nearest( const Point &pos ) const { return getNearestElement( pos ).second; }
double MeshGeometry::distance( const Point &pos, const Point &dir ) const
{
    // Update cached data if position moved
    if ( d_pos_hash != d_mesh->positionHash() )
        initializePosition();
    auto x = pos;
    // Update cached data if position moved
    if ( d_pos_hash != d_mesh->positionHash() )
        initializePosition();
    const auto type = d_mesh->getGeomType();
    while ( true ) {
        // Find the nearest node
        double d = 0;
        size_t i = d_tree.find_nearest( x.data(), &d );
        // For each parent element identify if we intersect the element
        for ( const auto &parent : d_mesh->getElementParents( d_nodes[i], type ) ) {
            // Find the distance to the element (if we intersect)
            double d = parent.distance( pos, dir );
            if ( d < 1e100 )
                return d;
        }
        // Nobody intersected, move the distance to the nearest node and repeat
        x += d * dir;
    }
    return 0;
}
bool MeshGeometry::inside( const Point &pos ) const
{
    double dist = distance( pos, pos - centroid() );
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
    auto elem = getNearestElement( x ).first;
    for ( auto id : d_surfaceIds ) {
        if ( elem.isInBlock( id ) )
            return id;
    }
    return 0;
}
Point MeshGeometry::surfaceNorm( const Point &x ) const
{
    auto elem = getNearestElement( x ).first;
    return elem.norm();
}


/********************************************************
 * Get the centroid/box                                  *
 ********************************************************/
Point MeshGeometry::centroid() const
{
    if ( d_pos_hash != d_mesh->positionHash() )
        initializePosition(); // Update cached data if position moved
    return d_centroid;
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
std::pair<AMP::Mesh::MeshElement, Point> MeshGeometry::getNearestElement( const Point &x ) const
{
    // Update cached data if position moved
    if ( d_pos_hash != d_mesh->positionHash() )
        initializePosition();
    // Get the nearest node
    size_t i = d_tree.find_nearest( x.data() );
    // Search each parent to determine the closest element
    double d = 1e200;
    Point p;
    AMP::Mesh::MeshElement elem;
    const auto type = d_mesh->getGeomType();
    for ( const auto &parent : d_mesh->getElementParents( d_nodes[i], type ) ) {
        auto p2 = parent.nearest( x );
        auto d2 = abs( p - p2 );
        if ( d2 < d ) {
            d    = d2;
            p    = p2;
            elem = parent;
        }
    }
    return std::make_pair( elem, p );
}


/********************************************************
 * Get the nearest element                               *
 ********************************************************/
bool MeshGeometry::isConvex() const
{
    bool is_convex = true;
    auto [lb, ub]  = box();
    double tol     = 1e-6 * abs( ub - lb );
    for ( const auto &elem : d_mesh->getIterator( d_mesh->getGeomType() ) ) {
        // Get the normal to the plane of the element
        auto n = elem.norm();
        // Get a point in the plane
        auto a = elem.getElements( AMP::Mesh::GeomType::Vertex )[0].coord();
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


} // namespace AMP::Geometry
