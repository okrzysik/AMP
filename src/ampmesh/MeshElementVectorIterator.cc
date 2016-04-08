#include "ampmesh/MeshElementVectorIterator.h"
#include "ampmesh/MeshElement.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Mesh {


// Create a unique id for this class
static unsigned int MultiVectorIteratorTypeID = TYPE_HASH( MultiVectorIterator );

/********************************************************
* Constructors                                          *
********************************************************/
MultiVectorIterator::MultiVectorIterator()
{
    typeID     = MultiVectorIteratorTypeID;
    iterator   = nullptr;
    d_elements = AMP::shared_ptr<std::vector<MeshElement>>();
    d_pos      = 0;
}
MultiVectorIterator::MultiVectorIterator( AMP::shared_ptr<std::vector<MeshElement>> elements,
                                          size_t pos )
{
    typeID     = MultiVectorIteratorTypeID;
    iterator   = nullptr;
    d_elements = elements;
    d_pos      = pos;
}
MultiVectorIterator::MultiVectorIterator( const std::vector<MeshElement> &elements, size_t pos )
{
    typeID   = MultiVectorIteratorTypeID;
    iterator = nullptr;
    d_elements.reset( new std::vector<MeshElement>( elements ) );
    d_pos = pos;
}
MultiVectorIterator::MultiVectorIterator( const MultiVectorIterator &rhs )
    : MeshIterator() // Note: we never want to call the base copy constructor
{
    typeID     = MultiVectorIteratorTypeID;
    iterator   = nullptr;
    d_elements = rhs.d_elements;
    d_pos      = rhs.d_pos;
}
MultiVectorIterator &MultiVectorIterator::operator=( const MultiVectorIterator &rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    this->typeID     = MultiVectorIteratorTypeID;
    this->iterator   = nullptr;
    this->d_elements = rhs.d_elements;
    this->d_pos      = rhs.d_pos;
    return *this;
}


/********************************************************
* Function to clone the iterator                        *
********************************************************/
MeshIterator *MultiVectorIterator::clone() const { return new MultiVectorIterator( *this ); }


/********************************************************
* De-constructor                                        *
********************************************************/
MultiVectorIterator::~MultiVectorIterator() {}


/********************************************************
* Return an iterator to the beginning or end            *
********************************************************/
MeshIterator MultiVectorIterator::begin() const { return MultiVectorIterator( d_elements, 0 ); }
MeshIterator MultiVectorIterator::end() const
{
    return MultiVectorIterator( d_elements, d_elements->size() );
}


/********************************************************
* Return the number of elements in the iterator         *
********************************************************/
size_t MultiVectorIterator::size() const
{
    if ( d_elements.get() == nullptr )
        return 0;
    return d_elements->size();
}
size_t MultiVectorIterator::position() const { return d_pos; }


/********************************************************
* Increment/Decrement the iterator                      *
********************************************************/
MeshIterator &MultiVectorIterator::operator++()
{
    // Prefix increment (increment and return this)
    d_pos++;
    return *this;
}
MeshIterator MultiVectorIterator::operator++( int )
{
    // Postfix increment (increment and return temporary object)
    MultiVectorIterator tmp( *this ); // Create a temporary variable
    this->operator++();               // apply operator
    return tmp;                       // return temporary result
}
MeshIterator &MultiVectorIterator::operator--()
{
    // Prefix decrement (increment and return this)
    d_pos--;
    return *this;
}
MeshIterator MultiVectorIterator::operator--( int )
{
    // Postfix decrement (increment and return temporary object)
    MultiVectorIterator tmp( *this ); // Create a temporary variable
    --( *this );                      // apply operator
    return tmp;                       // return temporary result
}


/********************************************************
* Random access incrementors                            *
********************************************************/
MeshIterator MultiVectorIterator::operator+( int n ) const
{
    MultiVectorIterator tmp( *this ); // Create a temporary iterator
    tmp.operator+=( n );              // Increment temporary iterator
    return tmp;
}
MeshIterator &MultiVectorIterator::operator+=( int n )
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
bool MultiVectorIterator::operator==( const MeshIterator &rhs ) const
{
    const MultiVectorIterator *rhs2 = nullptr;
    // Convert rhs to a MultiVectorIterator* so we can access the base class members
    const MultiVectorIterator *tmp = reinterpret_cast<const MultiVectorIterator *>( &rhs );
    if ( typeid( rhs ) == typeid( MultiVectorIterator ) ) {
        rhs2 = tmp; // We can safely cast rhs to a MultiVectorIterator
    } else if ( tmp->typeID == MultiVectorIteratorTypeID ) {
        rhs2 = tmp; // We can safely cast rhs.iterator to a MultiVectorIterator
    } else if ( ( reinterpret_cast<MultiVectorIterator *>( tmp->iterator ) )->typeID ==
                MultiVectorIteratorTypeID ) {
        rhs2 = reinterpret_cast<MultiVectorIterator *>( tmp->iterator );
    }
    // Perform direct comparisions if we are dealing with two MultiVectorIterators
    if ( rhs2 != nullptr ) {
        // Check that we are at the same position
        if ( d_pos != rhs2->d_pos )
            return false;
        // Check if we both arrays are the same memory address
        if ( d_elements.get() == rhs2->d_elements.get() )
            return true;
        // If we are dealing with different arrays, check that the are the same size and values
        if ( d_elements->size() != rhs2->d_elements->size() )
            return false;
        bool elements_match = true;
        for ( size_t i = 0; i < d_elements->size(); i++ ) {
            if ( d_elements->operator[]( i ) != rhs2->d_elements->operator[]( i ) )
                elements_match = false;
        }
        return elements_match;
    }
    /* We are comparing a MultiVectorIterator to an arbitrary iterator
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
    MeshIterator iterator = rhs.begin();
    bool elements_match   = true;
    for ( size_t i = 0; i < d_elements->size(); i++ ) {
        if ( iterator->globalID() != ( d_elements->operator[]( i ) ).globalID() )
            elements_match = false;
        ++iterator;
    }
    return elements_match;
}
bool MultiVectorIterator::operator!=( const MeshIterator &rhs ) const
{
    return !( ( *this ) == rhs );
}


/********************************************************
* Dereference the iterator to get the element           *
********************************************************/
MeshElement &MultiVectorIterator::operator*()
{
    if ( d_pos >= d_elements->size() )
        AMP_ERROR( "Invalid dereference (iterator is out of range" );
    return d_elements->operator[]( d_pos );
}
const MeshElement &MultiVectorIterator::operator*() const
{
    if ( d_pos >= d_elements->size() )
        AMP_ERROR( "Invalid dereference (iterator is out of range" );
    return d_elements->operator[]( d_pos );
}
MeshElement *MultiVectorIterator::operator->()
{
    if ( d_pos >= d_elements->size() )
        AMP_ERROR( "Invalid dereference (iterator is out of range" );
    return &( d_elements->operator[]( d_pos ) );
}
const MeshElement *MultiVectorIterator::operator->() const
{
    if ( d_pos >= d_elements->size() )
        AMP_ERROR( "Invalid dereference (iterator is out of range" );
    return &( d_elements->operator[]( d_pos ) );
}


} // Mesh namespace
} // AMP namespace
