#include "AMP/ampmesh/MeshElement.h"

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


} // namespace Mesh
} // namespace AMP
