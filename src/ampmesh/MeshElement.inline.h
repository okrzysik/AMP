#ifndef included_AMP_MeshElement_inline
#define included_AMP_MeshElement_inline

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
MeshElement::MeshElement( const MeshElement &rhs )
    : typeID( MeshElementTypeID ), element( nullptr ), d_globalID( rhs.d_globalID )
{
    if ( rhs.element == nullptr && rhs.typeID == MeshElementTypeID ) {
        element = nullptr;
    } else if ( rhs.typeID != MeshElementTypeID ) {
        element = rhs.clone();
    } else {
        element = rhs.element->clone();
    }
}
MeshElement::MeshElement( MeshElement &&rhs )
    : typeID( MeshElementTypeID ), element( rhs.element ), d_globalID( rhs.d_globalID )
{
    if ( rhs.typeID != MeshElementTypeID )
        element = rhs.clone();
    rhs.element = nullptr;
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
MeshElement &MeshElement::operator=( MeshElement &&rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    if ( element != nullptr ) {
        // Delete the existing element
        delete element;
        element = nullptr;
    }
    typeID     = MeshElementTypeID;
    d_globalID = rhs.d_globalID;
    std::swap( element, rhs.element );
    if ( rhs.typeID != MeshElementTypeID )
        element = rhs.clone();
    return *this;
}


/********************************************************
 * Destructor                                            *
 ********************************************************/
MeshElement::~MeshElement()
{
    if ( element != nullptr )
        delete element;
    element = nullptr;
}


/********************************************************
 * Is the element null                                   *
 ********************************************************/
bool MeshElement::isNull() const { return typeID == MeshElementTypeID && element == nullptr; }


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
inline MeshElement *MeshElement::getRawElement() { return element == nullptr ? this : element; }
inline const MeshElement *MeshElement::getRawElement() const
{
    return element == nullptr ? this : element;
}


/********************************************************
 * Function to check if a point is within an element     *
 ********************************************************/
bool MeshElement::containsPoint( const std::vector<double> &pos, double TOL ) const
{
    if ( element != nullptr )
        return element->containsPoint( pos, TOL );
    if ( d_globalID.type() == GeomType::Vertex ) {
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
 * Function to return basic info                         *
 ********************************************************/
inline const MeshElementID &MeshElement::globalID() const
{
    return element == nullptr ? d_globalID : element->d_globalID;
}
inline GeomType MeshElement::elementType() const { return globalID().type(); }
inline std::string MeshElement::elementClass() const
{
    return element == nullptr ? std::string( "MeshElement" ) : element->elementClass();
}


/********************************************************
 * Functions that are wrappers to an anvanced version    *
 ********************************************************/
inline std::vector<double> MeshElement::coord() const
{
    size_t N = 10;
    double x[10];
    ( element != nullptr ? element : this )->coord( N, x );
    return std::vector<double>( x, x + N );
}
inline double MeshElement::coord( int i ) const
{
    size_t N = 10;
    double x[10];
    ( element != nullptr ? element : this )->coord( N, x );
    return x[i];
}
inline std::vector<double> MeshElement::centroid() const
{
    size_t N = 10;
    double x[10];
    ( element != nullptr ? element : this )->centroid( N, x );
    return std::vector<double>( x, x + N );
}
inline std::vector<MeshElement> MeshElement::getElements( const GeomType type ) const
{
    std::vector<MeshElement> elements;
    ( element != nullptr ? element : this )->getElements( type, elements );
    return elements;
}
inline std::vector<MeshElement::shared_ptr> MeshElement::getNeighbors() const
{
    std::vector<MeshElement::shared_ptr> neighbors;
    ( element != nullptr ? element : this )->getNeighbors( neighbors );
    return neighbors;
}


/********************************************************
 * Functions that aren't implimented for the base class  *
 ********************************************************/
double MeshElement::volume() const
{
    if ( element == nullptr )
        AMP_ERROR( "volume is not implimented for the base class (" + elementClass() + ")" );
    return element->volume();
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


} // namespace Mesh
} // namespace AMP

#endif
