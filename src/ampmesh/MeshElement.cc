#include "AMP/ampmesh/MeshElement.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"

namespace AMP {
namespace Mesh {


/********************************************************
 * Function to return basic info         *
 ********************************************************/
std::string MeshElement::elementClass() const
{
    return element == nullptr ? std::string( "MeshElement" ) : element->elementClass();
}


/********************************************************
 * Function to return the centroid of an element         *
 ********************************************************/
Point MeshElement::centroid() const
{
    if ( element != nullptr )
        return element->centroid();
    if ( globalID().type() == GeomType::Vertex )
        return coord();
    std::vector<MeshElement> nodes;
    ( element != nullptr ? element : this )->getElements( GeomType::Vertex, nodes );
    AMP_ASSERT( !nodes.empty() );
    auto center = nodes[0].coord();
    for ( size_t i = 1; i < nodes.size(); i++ ) {
        auto pos = nodes[i].coord();
        for ( size_t j = 0; j < center.size(); j++ )
            center[j] += pos[j];
    }
    for ( size_t j = 0; j < center.size(); j++ )
        center[j] /= nodes.size();
    return center;
}


/********************************************************
 * Return the neighbors/elements                         *
 ********************************************************/
void MeshElement::getElements( const GeomType type, std::vector<MeshElement> &elements ) const
{
    if ( element == nullptr )
        AMP_ERROR( "getElements is not implimented for the base class (" + elementClass() + ")" );
    element->getElements( type, elements );
}
void MeshElement::getElementsID( const GeomType type, std::vector<MeshElementID> &ID ) const
{
    if ( element != nullptr )
        return element->getElementsID( type, ID );
    std::vector<MeshElement> elements;
    this->getElements( type, elements );
    ID.resize( elements.size() );
    for ( size_t i = 0; i < elements.size(); i++ )
        ID[i] = elements[i].globalID();
}
void MeshElement::getNeighbors( std::vector<MeshElement::shared_ptr> &neighbors ) const
{
    if ( element == nullptr )
        AMP_ERROR( "getNeighbors is not implimented for the base class (" + elementClass() + ")" );
    element->getNeighbors( neighbors );
}


/********************************************************
 * Function to check if a point is within an element     *
 ********************************************************/
bool MeshElement::containsPoint( const Point &pos, double TOL ) const
{
    if ( element != nullptr )
        return element->containsPoint( pos, TOL );
    if ( globalID().type() == GeomType::Vertex ) {
        // double dist = 0.0;
        auto point   = this->coord();
        double dist2 = 0.0;
        for ( size_t i = 0; i < point.size(); i++ )
            dist2 += ( point[i] - pos[i] ) * ( point[i] - pos[i] );
        return dist2 <= TOL * TOL;
    }
    AMP_ERROR( "containsPoint is not finished for default elements yet" );
    return false;
}


/********************************************************
 * Functions that aren't implimented for the base class  *
 ********************************************************/
Point MeshElement::coord() const
{
    if ( element == nullptr )
        AMP_ERROR( "coord is not implimented for the base class (" + elementClass() + ")" );
    return element->coord();
}
double MeshElement::volume() const
{
    if ( element == nullptr )
        AMP_ERROR( "volume is not implimented for the base class (" + elementClass() + ")" );
    return element->volume();
}
Point MeshElement::norm() const
{
    if ( element == nullptr )
        AMP_ERROR( "norm is not implimented for the base class (" + elementClass() + ")" );
    return element->norm();
}
Point MeshElement::nearest( const Point &pos ) const
{
    if ( element == nullptr )
        AMP_ERROR( "nearest is not implimented for the base class (" + elementClass() + ")" );
    return element->nearest( pos );
}
double MeshElement::distance( const Point &pos, const Point &dir ) const
{
    if ( element == nullptr )
        AMP_ERROR( "distance is not implimented for the base class (" + elementClass() + ")" );
    return element->distance( pos, dir );
}
bool MeshElement::isOnSurface() const
{
    if ( element == nullptr )
        AMP_ERROR( "isOnSurface is not implimented for the base class (" + elementClass() + ")" );
    return element->isOnSurface();
}
bool MeshElement::isOnBoundary( int id ) const
{
    if ( element == nullptr )
        AMP_ERROR( "isOnBoundary is not implimented for the base class (" + elementClass() + ")" );
    return element->isOnBoundary( id );
}
bool MeshElement::isInBlock( int id ) const
{
    if ( element == nullptr )
        AMP_ERROR( "isInBlock is not implimented for the base class (" + elementClass() + ")" );
    return element->isInBlock( id );
}
unsigned int MeshElement::globalOwnerRank() const
{
    if ( element == nullptr )
        AMP_ERROR( "globalOwnerRank is not implimented for the base class (" + elementClass() +
                   ")" );
    return element->globalOwnerRank();
}
MeshElementID MeshElement::globalID() const
{
    if ( element == nullptr )
        return MeshElementID();
    return element->globalID();
}


/********************************************************
 * Return points in the volume at the given resolution   *
 ********************************************************/
std::vector<Point> MeshElement::sample( double dx ) const
{
    auto type = globalID().type();
    if ( type == AMP::Mesh::GeomType::Vertex )
        return std::vector<Point>( 1, coord() );
    // Get the nodes
    auto nodes = getElements( AMP::Mesh::GeomType::Vertex );
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
        int n    = N;
        for ( int i = 0; i < n; i++ )
            p.push_back( x[0] + ( 0.5 * ( N - n ) + n ) * dx * ( x[1] - x[0] ) );
    } else if ( x.size() == 3u ) {
        // We are dealing with a triangle
        p = AMP::Geometry::GeometryHelpers::subdivide( { x[0], x[1], x[2] }, dx );
    } else if ( type == AMP::Mesh::GeomType::Face ) {
        // Get the normal
        auto n      = norm();
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
                if ( containsPoint( p1 ) )
                    p.push_back( p1 );
            }
        }
    } else {
        AMP_ERROR( "Not finished" );
    }
    return p;
}


} // namespace Mesh
} // namespace AMP
