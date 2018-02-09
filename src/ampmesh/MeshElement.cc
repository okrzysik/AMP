#include "AMP/ampmesh/MeshElement.h"

namespace AMP {
namespace Mesh {


/********************************************************
 * Function to return the centroid of an element         *
 ********************************************************/
void MeshElement::centroid( size_t &N, double *center ) const
{
    if ( element != nullptr )
        return element->centroid( N, center );
    if ( d_globalID.type() == GeomType::Vertex )
        return coord( N, center );
    std::vector<MeshElement> nodes;
    ( element != nullptr ? element : this )->getElements( GeomType::Vertex, nodes );
    AMP_ASSERT( !nodes.empty() );
    nodes[0].coord( N, center );
    for ( size_t i = 1; i < nodes.size(); i++ ) {
        double pos[10];
        nodes[i].coord( N, pos );
        for ( size_t j = 0; j < N; j++ )
            center[j] += pos[j];
    }
    for ( size_t j = 0; j < N; j++ )
        center[j] /= nodes.size();
}


/********************************************************
 * Function to return the coordinates of an element      *
 ********************************************************/
void MeshElement::coord( size_t &N, double *x ) const
{
    if ( element == nullptr )
        AMP_ERROR( "coord is not implimented for the base class (" + elementClass() + ")" );
    return element->coord( N, x );
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


} // namespace Mesh
} // namespace AMP
