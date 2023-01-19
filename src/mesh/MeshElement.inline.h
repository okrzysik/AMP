#ifndef included_AMP_MeshElement_inline
#define included_AMP_MeshElement_inline

#include "AMP/utils/UtilityMacros.h"

namespace AMP::Mesh {


/********************************************************
 * Constructors                                          *
 ********************************************************/
constexpr auto MeshElementHash = AMP::getTypeID<MeshElement>().hash;
static_assert( MeshElementHash != 0 );
MeshElement::MeshElement() : d_typeHash( MeshElementHash ), d_element( nullptr ) {}
MeshElement::MeshElement( const MeshElement &rhs )
    : d_typeHash( MeshElementHash ), d_element( nullptr )
{
    if ( rhs.d_element == nullptr && rhs.d_typeHash == MeshElementHash ) {
        d_element = nullptr;
    } else if ( rhs.d_typeHash != MeshElementHash ) {
        d_element = rhs.clone();
    } else {
        d_element = rhs.d_element->clone();
    }
}
MeshElement::MeshElement( MeshElement &&rhs )
    : d_typeHash( MeshElementHash ), d_element( rhs.d_element )
{
    if ( rhs.d_typeHash != MeshElementHash )
        d_element = rhs.clone();
    rhs.d_element = nullptr;
}
MeshElement &MeshElement::operator=( const MeshElement &rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    if ( d_element != nullptr ) {
        // Delete the existing d_element
        delete d_element;
        d_element = nullptr;
    }
    d_typeHash = MeshElementHash;
    if ( rhs.d_element == nullptr && rhs.d_typeHash == MeshElementHash ) {
        d_element = nullptr;
    } else if ( rhs.d_typeHash != MeshElementHash ) {
        d_element = rhs.clone();
    } else {
        d_element = rhs.d_element->clone();
    }
    return *this;
}
MeshElement &MeshElement::operator=( MeshElement &&rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    if ( d_element != nullptr ) {
        // Delete the existing d_element
        delete d_element;
        d_element = nullptr;
    }
    d_typeHash = MeshElementHash;
    std::swap( d_element, rhs.d_element );
    if ( rhs.d_typeHash != MeshElementHash )
        d_element = rhs.clone();
    return *this;
}
MeshElement::MeshElement( MeshElement *rhs ) : d_typeHash( MeshElementHash ), d_element( nullptr )
{
    if ( rhs->d_element ) {
        std::swap( d_element, rhs->d_element );
        delete rhs;
    } else {
        d_element = rhs;
    }
}


/********************************************************
 * Destructor                                            *
 ********************************************************/
MeshElement::~MeshElement()
{
    if ( d_element != nullptr )
        delete d_element;
    d_element = nullptr;
}


/********************************************************
 * Is the d_element null                                   *
 ********************************************************/
bool MeshElement::isNull() const { return d_typeHash == MeshElementHash && d_element == nullptr; }


/********************************************************
 * Function to clone the d_element                         *
 ********************************************************/
MeshElement *MeshElement::clone() const
{
    if ( d_element == nullptr )
        return new MeshElement();
    else
        AMP_ERROR( "clone must instantiated by the derived class" );
    return nullptr;
}


/********************************************************
 * Function to get the raw d_element                       *
 ********************************************************/
inline MeshElement *MeshElement::getRawElement() { return d_element == nullptr ? this : d_element; }
inline const MeshElement *MeshElement::getRawElement() const
{
    return d_element == nullptr ? this : d_element;
}


/********************************************************
 * Function to check if a point is within an d_element     *
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
    std::vector<MeshElement> d_elements;
    ( d_element != nullptr ? d_element : this )->getElements( type, d_elements );
    return d_elements;
}
inline std::vector<std::shared_ptr<MeshElement>> MeshElement::getNeighbors() const
{
    std::vector<std::shared_ptr<MeshElement>> neighbors;
    ( d_element != nullptr ? d_element : this )->getNeighbors( neighbors );
    return neighbors;
}


/****************************************************************
 * Stream operator                                               *
 ****************************************************************/
std::ostream &operator<<( std::ostream &out, const AMP::Mesh::MeshElement &x );


} // namespace AMP::Mesh


#endif
