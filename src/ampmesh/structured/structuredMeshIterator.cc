#include "ampmesh/structured/structuredMeshIterator.h"
#include "ampmesh/structured/structuredMeshElement.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Mesh {


// Create a unique id for this class
static unsigned int structuredMeshIteratorTypeID = TYPE_HASH( structuredMeshIterator );

// unused global variable to prevent compiler warning
static MeshElement nullElement;


/********************************************************
* Constructors                                          *
********************************************************/
structuredMeshIterator::structuredMeshIterator()
{
    typeID   = structuredMeshIteratorTypeID;
    iterator = nullptr;
    d_pos    = 0;
    d_mesh   = nullptr;
}
structuredMeshIterator::structuredMeshIterator(
    AMP::shared_ptr<std::vector<BoxMesh::MeshElementIndex>> elements,
    const AMP::Mesh::BoxMesh *mesh,
    size_t pos )
{
    typeID     = structuredMeshIteratorTypeID;
    iterator   = nullptr;
    d_pos      = pos;
    d_mesh     = mesh;
    d_elements = elements;
}
structuredMeshIterator::structuredMeshIterator( const structuredMeshIterator &rhs ) : MeshIterator()
{
    typeID     = structuredMeshIteratorTypeID;
    iterator   = nullptr;
    d_pos      = rhs.d_pos;
    d_mesh     = rhs.d_mesh;
    d_elements = rhs.d_elements;
}
structuredMeshIterator &structuredMeshIterator::operator=( const structuredMeshIterator &rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    this->typeID     = structuredMeshIteratorTypeID;
    this->iterator   = nullptr;
    this->d_pos      = rhs.d_pos;
    this->d_mesh     = rhs.d_mesh;
    this->d_elements = rhs.d_elements;
    return *this;
}


/********************************************************
* Function to clone the iterator                        *
********************************************************/
MeshIterator *structuredMeshIterator::clone() const { return new structuredMeshIterator( *this ); }


/********************************************************
* De-constructor                                        *
********************************************************/
structuredMeshIterator::~structuredMeshIterator() {}


/********************************************************
* Return an iterator to the beginning or end            *
********************************************************/
MeshIterator structuredMeshIterator::begin() const
{
    return structuredMeshIterator( d_elements, d_mesh, 0 );
}
MeshIterator structuredMeshIterator::end() const
{
    return structuredMeshIterator( d_elements, d_mesh, d_elements->size() );
}


/********************************************************
* Return the number of elements in the iterator         *
********************************************************/
size_t structuredMeshIterator::size() const { return d_elements->size(); }
size_t structuredMeshIterator::position() const { return d_pos; }


/********************************************************
* Increment/Decrement the iterator                      *
********************************************************/
MeshIterator &structuredMeshIterator::operator++()
{
    // Prefix increment (increment and return this)
    d_pos++;
    if ( d_pos > d_elements->size() )
        d_pos = d_elements->size();
    return *this;
}
MeshIterator structuredMeshIterator::operator++( int )
{
    // Postfix increment (increment and return temporary object)
    structuredMeshIterator tmp( *this ); // Create a temporary variable
    this->operator++();                  // apply operator
    return tmp;                          // return temporary result
}
MeshIterator &structuredMeshIterator::operator--()
{
    // Prefix decrement (increment and return this)
    if ( d_pos != 0 )
        d_pos--;
    return *this;
}
MeshIterator structuredMeshIterator::operator--( int )
{
    // Postfix decrement (increment and return temporary object)
    structuredMeshIterator tmp( *this ); // Create a temporary variable
    --( *this );                         // apply operator
    return tmp;                          // return temporary result
}


/********************************************************
* Random access incrementors                            *
********************************************************/
MeshIterator structuredMeshIterator::operator+( int n ) const
{
    structuredMeshIterator tmp( *this ); // Create a temporary iterator
    tmp.operator+=( n );                 // Increment temporary iterator
    return tmp;
}
MeshIterator &structuredMeshIterator::operator+=( int n )
{
    if ( n >= 0 ) { // increment *this
        size_t n2 = static_cast<size_t>( n );
        if ( d_pos + n2 > d_elements->size() )
            AMP_ERROR( "Iterated past end of iterator" );
        d_pos += n2;
    } else { // decrement *this
        size_t n2 = static_cast<size_t>( -n );
        if ( n2 > d_pos )
            AMP_ERROR( "Iterated past beginning of iterator" );
        d_pos -= n2;
    }
    return *this;
}


/********************************************************
* Compare two iterators                                 *
********************************************************/
bool structuredMeshIterator::operator==( const MeshIterator &rhs ) const
{
    structuredMeshIterator *rhs2 = nullptr;
    structuredMeshIterator *tmp =
        (structuredMeshIterator *) &rhs; // Convert rhs to a structuredMeshIterator* so we can
                                         // access the base class members
    if ( typeid( rhs ) == typeid( structuredMeshIterator ) ) {
        rhs2 = tmp; // We can safely cast rhs to a structuredMeshIterator
    } else if ( tmp->typeID == structuredMeshIteratorTypeID ) {
        rhs2 = tmp; // We can safely cast rhs.iterator to a structuredMeshIterator
    } else if ( ( (structuredMeshIterator *) tmp->iterator )->typeID ==
                structuredMeshIteratorTypeID ) {
        rhs2 = (structuredMeshIterator *) tmp->iterator;
    }
    // Perform direct comparisions if we are dealing with two structuredMeshIterators
    if ( rhs2 != nullptr ) {
        if ( d_pos != rhs2->d_pos )
            return false;
        if ( d_elements->size() != rhs2->d_elements->size() )
            return false;
        if ( d_elements.get() == d_elements.get() )
            return true;
        for ( size_t i = 0; i < d_elements->size(); i++ ) {
            if ( d_elements->operator[]( i ) != rhs2->d_elements->operator[]( i ) )
                return false;
        }
        return true;
    }
    /* We are comparing a structuredMeshIterator to an arbitrary iterator
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
    AMP::Mesh::MeshIterator iterator = rhs.begin();
    for ( size_t i = 0; i < d_elements->size(); i++ ) {
        structuredMeshElement *elem =
            dynamic_cast<structuredMeshElement *>( iterator->getRawElement() );
        if ( elem == nullptr )
            return false;
        if ( elem->d_index != d_elements->operator[]( i ) )
            return false;
        ++iterator;
    }
    return true;
}
bool structuredMeshIterator::operator!=( const MeshIterator &rhs ) const
{
    return !( ( *this ) == rhs );
}


/********************************************************
* Dereference the iterator to get the element           *
********************************************************/
MeshElement &structuredMeshIterator::operator*()
{
    this->operator->(); // Initialize d_cur_element
    return d_cur_element;
}
MeshElement *structuredMeshIterator::operator->()
{
    d_cur_element = structuredMeshElement( d_elements->operator[]( d_pos ), d_mesh );
    return &d_cur_element;
}


} // Mesh namespace
} // AMP namespace
