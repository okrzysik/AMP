#include "AMP/ampmesh/MeshUtilities.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/utils/NearestPairSearch.h"
#include "AMP/utils/kdtree.h"

namespace AMP::Mesh {


// Helper function to check that no two points are too close
static inline double nearest( const std::vector<Point> &x )
{
    if ( x.empty() )
        return 0;
    auto index = find_min_dist( x );
    auto dx    = x[index.first] - x[index.second];
    return dx.abs();
}


static inline double dist( const kdtree &tree, const Point &p )
{
    double d;
    tree.find_nearest( p.data(), &d );
    return d;
}


/********************************************************
 *  Get points in the mesh                               *
 ********************************************************/
std::tuple<std::vector<Point>, std::vector<MeshElementID>> sample( const Mesh &mesh, double dx )
{
    std::vector<Point> points;
    std::vector<MeshElementID> ids;
    const auto type = mesh.getGeomType();
    // Start by adding all the nodes to preserve the geometry
    for ( auto node : mesh.getIterator( AMP::Mesh::GeomType::Vertex ) ) {
        points.push_back( node.coord() );
        ids.push_back( node.globalID() );
    }
    double d0 = nearest( points );
    // Check if we are auto determining the resolution
    if ( dx == 0.0 )
        dx = d0;
    // Create a kdtree to quickly calculate distances
    kdtree tree( points );
    // Function to add points
    auto addPoints = [&points, &ids, &tree, dx]( const MeshIterator &it ) {
        for ( const auto &elem : it ) {
            auto points2 = sample( elem, dx );
            for ( const auto &p : points2 ) {
                if ( dist( tree, p ) > 0.8 * dx ) {
                    points.push_back( p );
                    ids.push_back( elem.globalID() );
                    tree.add( p.data() );
                }
            }
        }
    };
    // Add points along the edges
    if ( type >= GeomType::Edge )
        addPoints( mesh.getIterator( GeomType::Edge ) );
    // Add points along the faces
    if ( type >= GeomType::Face )
        addPoints( mesh.getIterator( GeomType::Face ) );
    // Add points along the volume
    if ( type >= GeomType::Volume )
        addPoints( mesh.getIterator( GeomType::Volume ) );
    return std::tie( points, ids );
}


/********************************************************
 * Return points in the element at the given resolution  *
 ********************************************************/
std::vector<Point> sample( const MeshElement &elem, double dx )
{
    auto type = elem.globalID().type();
    if ( type == AMP::Mesh::GeomType::Vertex )
        return std::vector<Point>( 1, elem.coord() );
    // Get the nodes
    auto nodes = elem.getElements( AMP::Mesh::GeomType::Vertex );
    std::vector<Point> x( nodes.size() );
    for ( size_t i = 0; i < nodes.size(); i++ )
        x[i] = nodes[i].coord();
    // Check if we are dealing with a volume (in the coordinate space)
    if ( static_cast<int>( type ) == x[0].ndim() ) {
        // Create a uniform grid
        AMP_ERROR( "Not finished" );
    }
    // Code for the different object types
    std::vector<Point> p;
    if ( type == AMP::Mesh::GeomType::Edge ) {
        // We are dealing with an edge (easy)
        AMP_ASSERT( x.size() == 2u );
        double d = ( x[1] - x[0] ).abs();
        double N = d / dx;
        int n    = N + 1;
        auto dx2 = 1.0 / n * ( x[1] - x[0] );
        for ( int i = 0; i <= n; i++ )
            p.push_back( x[0] + static_cast<double>( i ) * dx2 );
    } else if ( x.size() == 3u ) {
        // We are dealing with a triangle
        p = AMP::Geometry::GeometryHelpers::subdivide( { x[0], x[1], x[2] }, dx );
    } else if ( type == AMP::Mesh::GeomType::Face ) {
        // Get the normal
        auto n      = elem.norm();
        auto center = x[0];
        for ( size_t i = 1; i < nodes.size(); i++ )
            center += x[i];
        center *= 1.0 / nodes.size();
        // Get two perpendicular unit vectors in the plane
        // Choose the furthest point for the first vector
        Point v1;
        double tmp = 0;
        for ( size_t i = 0; i < x.size(); i++ ) {
            auto v = x[i] - center;
            auto d = v.norm();
            if ( d > tmp ) {
                v1  = normalize( v );
                tmp = d;
            }
        }
        // Compute the second vector using the norm and the first
        auto v2 = normalize( cross( n, v1 ) );
        // Get the range in the unit vectors
        double range[4] = { 0, 0, 0, 0 };
        for ( size_t i = 0; i < x.size(); i++ ) {
            double d1 = dot( x[i] - center, v1 );
            double d2 = dot( x[i] - center, v2 );
            range[0]  = std::min( range[0], d1 );
            range[1]  = std::max( range[1], d1 );
            range[2]  = std::min( range[2], d2 );
            range[3]  = std::max( range[3], d2 );
        }
        // Create the points
        auto N1 = ( range[1] - range[0] ) / dx;
        auto N2 = ( range[3] - range[2] ) / dx;
        auto x1 = range[0] + 0.5 * ( N1 - floor( N1 ) );
        auto x2 = range[2] + 0.5 * ( N2 - floor( N2 ) );
        for ( double d1 = x1; d1 <= range[1]; d1 += dx ) {
            for ( double d2 = x2; d2 <= range[3]; d2 += dx ) {
                Point p1 = center + d1 * v1 + d2 * v2;
                if ( elem.containsPoint( p1 ) )
                    p.push_back( p1 );
            }
        }
    } else {
        AMP_ERROR( "Not finished" );
    }
    return p;
}

} // namespace AMP::Mesh
