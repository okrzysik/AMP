#include "AMP/ampmesh/triangle/TriangleMeshIterator.h"
#include "AMP/ampmesh/triangle/TriangleMesh.h"
#include "AMP/ampmesh/triangle/TriangleMeshElement.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Mesh {


// unused global variable to prevent compiler warning
static MeshElement nullElement;


/********************************************************
 * Create a unique id for each class                     *
 ********************************************************/
template<>
constexpr uint32_t TriangleMeshIterator<1, 1>::getTypeID()
{
    return AMP::Utilities::hash_char( "TriangleMeshIterator<1,1>" );
}
template<>
constexpr uint32_t TriangleMeshIterator<1, 2>::getTypeID()
{
    return AMP::Utilities::hash_char( "TriangleMeshIterator<1,2>" );
}
template<>
constexpr uint32_t TriangleMeshIterator<1, 3>::getTypeID()
{
    return AMP::Utilities::hash_char( "TriangleMeshIterator<1,3>" );
}
template<>
constexpr uint32_t TriangleMeshIterator<2, 2>::getTypeID()
{
    return AMP::Utilities::hash_char( "TriangleMeshIterator<2,2>" );
}
template<>
constexpr uint32_t TriangleMeshIterator<2, 3>::getTypeID()
{
    return AMP::Utilities::hash_char( "TriangleMeshIterator<2,3>" );
}
template<>
constexpr uint32_t TriangleMeshIterator<3, 3>::getTypeID()
{
    return AMP::Utilities::hash_char( "TriangleMeshIterator<3,3>" );
}


/********************************************************
 * Constructors                                          *
 ********************************************************/
template<size_t NG, size_t NP>
TriangleMeshIterator<NG, NP>::TriangleMeshIterator()
{
    d_typeID   = getTypeID();
    d_iterator = nullptr;
    d_size     = 0;
    d_pos      = -1;
    d_element  = &d_cur_element;
    d_mesh     = nullptr;
}
template<size_t NG, size_t NP>
TriangleMeshIterator<NG, NP>::TriangleMeshIterator(
    const AMP::Mesh::TriangleMesh<NG, NP> *mesh,
    AMP::shared_ptr<const std::vector<ElementID>> list,
    size_t pos )
{
    d_typeID   = getTypeID();
    d_iterator = nullptr;
    d_size     = 0;
    d_pos      = pos;
    d_element  = &d_cur_element;
    d_mesh     = mesh;
    d_list     = list;
    if ( list )
        d_size = list->size();
    d_cur_element =
        TriangleMeshElement<NG, NP>( MeshElementID( mesh->meshID(), ElementID() ), mesh );
    if ( d_pos < d_size )
        d_cur_element.resetElemId( d_list->operator[]( d_pos ) );
}
template<size_t NG, size_t NP>
TriangleMeshIterator<NG, NP>::TriangleMeshIterator( const TriangleMeshIterator &rhs )
    : MeshIterator()
{
    d_typeID      = rhs.d_typeID;
    d_iterator    = nullptr;
    d_size        = rhs.d_size;
    d_pos         = rhs.d_pos;
    d_mesh        = rhs.d_mesh;
    d_list        = rhs.d_list;
    d_element     = &d_cur_element;
    d_cur_element = rhs.d_cur_element;
}
template<size_t NG, size_t NP>
TriangleMeshIterator<NG, NP> &TriangleMeshIterator<NG, NP>::
operator=( const TriangleMeshIterator &rhs )
{
    if ( this == &rhs )
        return *this;
    d_typeID      = rhs.d_typeID;
    d_iterator    = nullptr;
    d_size        = rhs.d_size;
    d_pos         = rhs.d_pos;
    d_mesh        = rhs.d_mesh;
    d_list        = rhs.d_list;
    d_element     = &d_cur_element;
    d_cur_element = rhs.d_cur_element;
    return *this;
}


/********************************************************
 * Function to clone the iterator                        *
 ********************************************************/
template<size_t NG, size_t NP>
MeshIterator *TriangleMeshIterator<NG, NP>::clone() const
{
    return new TriangleMeshIterator( *this );
}


/********************************************************
 * Return an iterator to the beginning or end            *
 ********************************************************/
template<size_t NG, size_t NP>
MeshIterator TriangleMeshIterator<NG, NP>::begin() const
{
    return TriangleMeshIterator( d_mesh, d_list, 0 );
}
template<size_t NG, size_t NP>
MeshIterator TriangleMeshIterator<NG, NP>::end() const
{
    return TriangleMeshIterator( d_mesh, d_list, d_size );
}


/********************************************************
 * Increment/Decrement the iterator                      *
 ********************************************************/
template<size_t NG, size_t NP>
MeshIterator &TriangleMeshIterator<NG, NP>::operator++()
{
    // Prefix increment (increment and return this)
    d_pos++;
    if ( d_pos < d_size )
        d_cur_element.resetElemId( d_list->operator[]( d_pos ) );
    return *this;
}
template<size_t NG, size_t NP>
MeshIterator TriangleMeshIterator<NG, NP>::operator++( int )
{
    // Postfix increment (increment and return temporary object)
    TriangleMeshIterator tmp( *this ); // Create a temporary variable
    this->operator++();                // apply operator
    return tmp;                        // return temporary result
}
template<size_t NG, size_t NP>
MeshIterator &TriangleMeshIterator<NG, NP>::operator--()
{
    // Prefix decrement (increment and return this)
    AMP_INSIST( d_pos > 0, "Decrementing iterator past 0" );
    d_pos--;
    d_cur_element.resetElemId( d_list->operator[]( d_pos ) );
    return *this;
}
template<size_t NG, size_t NP>
MeshIterator TriangleMeshIterator<NG, NP>::operator--( int )
{
    // Postfix decrement (increment and return temporary object)
    TriangleMeshIterator tmp( *this ); // Create a temporary variable
    --( *this );                       // apply operator
    return tmp;                        // return temporary result
}


/********************************************************
 * Random access incrementors                            *
 ********************************************************/
template<size_t NG, size_t NP>
MeshIterator TriangleMeshIterator<NG, NP>::operator+( int n ) const
{
    TriangleMeshIterator tmp( *this ); // Create a temporary iterator
    tmp.operator+=( n );               // Increment temporary iterator
    return tmp;
}
template<size_t NG, size_t NP>
MeshIterator &TriangleMeshIterator<NG, NP>::operator+=( int n )
{
    // Check the input
    if ( n >= 0 ) {
        AMP_INSIST( d_pos + n <= d_size, "Iterated past end of iterator" );
    } else { // decrement *this
        AMP_INSIST( -n <= (int64_t) d_pos, "Iterated past beginning of iterator" );
    }
    // Perform the increment and return
    d_pos += n;
    if ( d_pos < d_size )
        d_cur_element.resetElemId( d_list->operator[]( d_pos ) );
    return *this;
}


/********************************************************
 * Compare two iterators                                 *
 ********************************************************/
template<size_t NG, size_t NP>
bool TriangleMeshIterator<NG, NP>::operator==( const MeshIterator &rhs ) const
{
    const TriangleMeshIterator *rhs2 = nullptr;
    // Convert rhs to a TriangleMeshIterator* so we can access the base class members
    auto *tmp = reinterpret_cast<const TriangleMeshIterator *>( &rhs );
    if ( tmp->d_typeID == getTypeID() ) {
        rhs2 = tmp; // We can safely cast rhs to a TriangleMeshIterator
    } else if ( tmp->d_iterator != nullptr ) {
        tmp = reinterpret_cast<const TriangleMeshIterator *>( tmp->d_iterator );
        if ( tmp->d_typeID == getTypeID() )
            rhs2 = tmp; // We can safely cast rhs.iterator to a TriangleMeshIterator
    }
    // Perform direct comparisions if we are dealing with two TriangleMeshIterators
    if ( rhs2 != nullptr )
        return *d_list == *rhs2->d_list;
    /* We are comparing a TriangleMeshIterator to an arbitrary iterator
     * The iterators are the same if they point to the same position and iterate
     * over the same elements in the same order
     */
    // Check the size
    if ( this->size() != rhs.size() )
        return false;
    // Check the current position
    if ( this->position() != rhs.position() )
        return false;
    // Check that the elements match
    MeshIterator it1    = this->begin();
    MeshIterator it2    = rhs.begin();
    bool elements_match = true;
    for ( size_t i = 0; i < it1.size(); ++i, ++it1, ++it2 ) {
        if ( it1->globalID() != it2->globalID() )
            elements_match = false;
    }
    return elements_match;
}
template<size_t NG, size_t NP>
bool TriangleMeshIterator<NG, NP>::operator!=( const MeshIterator &rhs ) const
{
    return !( ( *this ) == rhs );
}


/********************************************************
 *  Explicit instantiations of TriangleMeshIterator      *
 ********************************************************/
template class TriangleMeshIterator<1, 1>;
template class TriangleMeshIterator<1, 2>;
template class TriangleMeshIterator<1, 3>;
template class TriangleMeshIterator<2, 2>;
template class TriangleMeshIterator<2, 3>;
template class TriangleMeshIterator<3, 3>;


} // namespace Mesh
} // namespace AMP
