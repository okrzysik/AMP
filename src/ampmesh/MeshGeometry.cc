#include "AMP/ampmesh/MeshGeometry.h"
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshElement.h"
#include "AMP/ampmesh/MeshUtilities.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
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
    std::vector<AMP::Mesh::MeshPoint<double>> points;
    std::tie( points, d_ids ) = AMP::Mesh::sample( *d_mesh, 0 );
    d_tree                    = AMP::kdtree( points );
}

/********************************************************
 * Get the distance to the surface                       *
 ********************************************************/
Point MeshGeometry::nearest( const Point &pos ) const { return getNearestPoint( pos ).second; }
double MeshGeometry::distance( const Point &pos, const Point &dir ) const
{
    NULL_USE( pos );
    NULL_USE( dir );
    /*    // Update cached data if position moved
        if ( d_pos_hash != d_mesh->positionHash() )
            initializePosition();
        // Get the bounding box for the mesh
        std::array<double, 6> box = { 0, 0, 0, 0, 0, 0 };
        const auto box2           = d_mesh->getBoundingBox();
        for ( size_t i = 0; i < box2.size(); i++ )
            box[i] = box2[i];
        // Compute the distance
        auto x               = pos;
        const auto type      = d_mesh->getGeomType();
        constexpr double inf = std::numeric_limits<double>::infinity();
        while ( true ) {
            // Check that the updated ray will intersect the bounding box
            if ( AMP::Geometry::GeometryHelpers::distanceToBox( x, dir, box ) == inf )
                return inf;
            // Find the nearest node
            double d = 0;
            size_t i = d_tree.find_nearest( x.data(), &d );
            AMP_ASSERT( i < d_nodes.size() );
            // For each parent element identify if we intersect the element
            for ( const auto &parent : d_mesh->getElementParents( d_nodes[i], type ) ) {
                // Find the distance to the element (if we intersect)
                double d = parent.distance( pos, dir );
                if ( d < 1e100 )
                    return d;
            }
            // Nobody intersected, move the distance to the nearest node and repeat
            x += d * dir;
        }*/
    AMP_ERROR( "Not finished" );
    return 0;
}
bool MeshGeometry::inside( const Point &pos ) const
{
    // Update cached data if position moved
    if ( d_pos_hash != d_mesh->positionHash() )
        initializePosition();
    // Get the nearest elements
    auto elems = getNearestElements( pos );
    // If the nearest element is a surface, add the neighbors
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

    return test;
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
    auto elem = getNearestPoint( x ).first;
    for ( auto id : d_surfaceIds ) {
        if ( elem.isInBlock( id ) )
            return id;
    }
    return 0;
}
Point MeshGeometry::surfaceNorm( const Point &x ) const
{
    auto elem = getNearestPoint( x ).first;
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
std::vector<AMP::Mesh::MeshElement> MeshGeometry::getNearestElements( const Point &x ) const
{
    // Update cached data if position moved
    if ( d_pos_hash != d_mesh->positionHash() )
        initializePosition();
    // Get the nearest node/edge/face
    size_t i = d_tree.find_nearest( x.data() );
    auto id0 = d_ids[i];
    // Find the nearest point and parent element
    auto type = d_mesh->getGeomType();
    auto elem = d_mesh->getElement( id0 );
    if ( id0.type() == type )
        return { elem };
    else
        return d_mesh->getElementParents( elem, type );
}
std::pair<AMP::Mesh::MeshElement, Point> MeshGeometry::getNearestPoint( const Point &x ) const
{
    // Update cached data if position moved
    if ( d_pos_hash != d_mesh->positionHash() )
        initializePosition();
    // Get the nearest node/edge/face
    size_t i = d_tree.find_nearest( x.data() );
    auto id0 = d_ids[i];
    // Find the nearest point and parent element
    auto type = d_mesh->getGeomType();
    auto elem = d_mesh->getElement( id0 );
    if ( id0.type() == type ) {
        // The element is the parent type
        auto p = elem.nearest( x );
        return std::make_pair( elem, p );
    } else {
        // Check each parent to find the closest point
        double d = 1e200;
        Point p;
        AMP::Mesh::MeshElement elem2;
        auto parents = d_mesh->getElementParents( elem, type );
        AMP_INSIST( !parents.empty(), "Error with mesh, node has no parents" );
        for ( const auto &parent : parents ) {
            auto p2 = parent.nearest( x );
            auto d2 = ( x - p2 ).norm();
            if ( d2 < d ) {
                d     = d2;
                p     = p2;
                elem2 = parent;
            }
        }
        return std::make_pair( elem2, p );
    }
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


} // namespace AMP::Geometry
