#ifndef included_AMP_MeshElement_inline
#define included_AMP_MeshElement_inline

#include "AMP/utils/UtilityMacros.h"

namespace AMP::Mesh {


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
inline std::vector<std::unique_ptr<MeshElement>> MeshElement::getNeighbors() const
{
    std::vector<std::unique_ptr<MeshElement>> neighbors;
    ( d_element != nullptr ? d_element : this )->getNeighbors( neighbors );
    return neighbors;
}


/****************************************************************
 * Stream operator                                               *
 ****************************************************************/
std::ostream &operator<<( std::ostream &out, const AMP::Mesh::MeshElement &x );


} // namespace AMP::Mesh


#endif
