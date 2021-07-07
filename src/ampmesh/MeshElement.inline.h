#ifndef included_AMP_MeshElement_inline
#define included_AMP_MeshElement_inline

#include "AMP/utils/UtilityMacros.h"

namespace AMP {
namespace Mesh {


/********************************************************
 * Constructors                                          *
 ********************************************************/
MeshElement::MeshElement()
{
    typeID  = MeshElementTypeID;
    element = nullptr;
}
MeshElement::MeshElement( const MeshElement &rhs ) : typeID( MeshElementTypeID ), element( nullptr )
{
    if ( rhs.element == nullptr && rhs.typeID == MeshElementTypeID ) {
        element = nullptr;
    } else if ( rhs.typeID != MeshElementTypeID ) {
        element = rhs.clone();
    } else {
        element = rhs.element->clone();
    }
}
MeshElement::MeshElement( MeshElement &&rhs ) : typeID( MeshElementTypeID ), element( rhs.element )
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
    typeID = MeshElementTypeID;
    std::swap( element, rhs.element );
    if ( rhs.typeID != MeshElementTypeID )
        element = rhs.clone();
    return *this;
}
MeshElement::MeshElement( MeshElement *rhs ) : typeID( MeshElementTypeID ), element( nullptr )
{
    if ( rhs->element ) {
        std::swap( element, rhs->element );
        delete rhs;
    } else {
        element = rhs;
    }
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
inline bool MeshElement::containsPoint( const std::vector<double> &pos, double TOL ) const
{
    return containsPoint( Point( pos.size(), pos.data() ), TOL );
}


/********************************************************
 * Functions that are wrappers to an advanced version    *
 ********************************************************/
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


/****************************************************************
 * Stream operator                                               *
 ****************************************************************/
std::ostream &operator<<( std::ostream &out, const AMP::Mesh::MeshElement &x );


} // namespace Mesh
} // namespace AMP


#endif
