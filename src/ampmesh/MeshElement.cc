#include "ampmesh/MeshElement.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Mesh {


// Create a unique id for this class
static unsigned int MeshElementTypeID = TYPE_HASH( MeshElement );


/********************************************************
* Constructors                                          *
********************************************************/
MeshElement::MeshElement() : d_globalID()
{
    typeID  = MeshElementTypeID;
    element = nullptr;
}
MeshElement::MeshElement( const MeshElement &rhs ) : d_globalID()
{
    typeID  = MeshElementTypeID;
    element = nullptr;
    if ( rhs.element == nullptr && rhs.typeID == MeshElementTypeID ) {
        element = nullptr;
    } else if ( rhs.typeID != MeshElementTypeID ) {
        element = rhs.clone();
    } else {
        element = rhs.element->clone();
    }
    d_globalID = rhs.d_globalID;
}
MeshElement &MeshElement::operator=( const MeshElement &rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    if ( element != nullptr ) {
        // Delete the existing element
        delete element;
        element = nullptr;
    }
    typeID = MeshElementTypeID;
    if ( rhs.element == nullptr && rhs.typeID == MeshElementTypeID ) {
        element = nullptr;
    } else if ( rhs.typeID != MeshElementTypeID ) {
        element = rhs.clone();
    } else {
        element = rhs.element->clone();
    }
    d_globalID = rhs.d_globalID;
    return *this;
}


/********************************************************
* De-constructor                                        *
********************************************************/
MeshElement::~MeshElement()
{
    if ( element != nullptr )
        delete element;
    element = nullptr;
}


/********************************************************
* Function to clone the element                         *
********************************************************/
MeshElement *MeshElement::clone() const
{
    if ( element == nullptr )
        return new MeshElement();
    else
        AMP_ERROR( "clone must instantiated by the derived class" );
    return nullptr;
}


/********************************************************
* Function to get the raw element                       *
********************************************************/
MeshElement *MeshElement::getRawElement()
{
    if ( element == nullptr )
        return this;
    return element->getRawElement();
}
const MeshElement *MeshElement::getRawElement() const
{
    if ( element == nullptr )
        return this;
    return element->getRawElement();
}


/********************************************************
* Function to return the centroid of an element         *
********************************************************/
std::vector<double> MeshElement::centroid() const
{
    if ( element != nullptr )
        return element->centroid();
    if ( d_globalID.type() == Vertex )
        return coord();
    std::vector<MeshElement> nodes = getElements( Vertex );
    AMP_ASSERT( !nodes.empty() );
    std::vector<double> center = nodes[0].coord();
    for ( size_t i = 1; i < nodes.size(); i++ ) {
        std::vector<double> coord = nodes[i].coord();
        for ( size_t j = 0; j < center.size(); j++ )
            center[j] += coord[j];
    }
    for ( auto &elem : center )
        elem /= nodes.size();
    return center;
}


/********************************************************
* Function to check if a point is within an element     *
********************************************************/
bool MeshElement::containsPoint( const std::vector<double> &pos, double TOL ) const
{
    if ( element != nullptr )
        return element->containsPoint( pos, TOL );
    if ( d_globalID.type() == Vertex ) {
        // double dist = 0.0;
        std::vector<double> point = this->coord();
        double dist2              = 0.0;
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
std::vector<MeshElement> MeshElement::getElements( const GeomType type ) const
{
    if ( element == nullptr )
        AMP_ERROR( "getElements is not implimented for the base class" );
    return element->getElements( type );
}
std::vector<MeshElement::shared_ptr> MeshElement::getNeighbors() const
{
    if ( element == nullptr )
        AMP_ERROR( "getNeighbors is not implimented for the base class" );
    return element->getNeighbors();
}
double MeshElement::volume() const
{
    if ( element == nullptr )
        AMP_ERROR( "volume is not implimented for the base class" );
    return element->volume();
}
std::vector<double> MeshElement::coord() const
{
    if ( element == nullptr )
        AMP_ERROR( "coord is not implimented for the base class" );
    return element->coord();
}
double MeshElement::coord( int i ) const
{
    if ( element == nullptr )
        return this->coord()[i];
    return element->coord( i );
}
bool MeshElement::isOnSurface() const
{
    if ( element == nullptr )
        AMP_ERROR( "isOnSurface is not implimented for the base class" );
    return element->isOnSurface();
}
bool MeshElement::isOnBoundary( int id ) const
{
    if ( element == nullptr )
        AMP_ERROR( "isOnBoundary is not implimented for the base class" );
    return element->isOnBoundary( id );
}
bool MeshElement::isInBlock( int id ) const
{
    if ( element == nullptr )
        AMP_ERROR( "isInBlock is not implimented for the base class" );
    return element->isInBlock( id );
}


} // Mesh namespace
} // AMP namespace
