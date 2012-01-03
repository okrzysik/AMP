#include "ampmesh/MeshElementVectorIterator.h"
#include "ampmesh/MeshElement.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Mesh {


// Create a unique id for this class
static unsigned int MultiVectorIteratorTypeID = TYPE_HASH(MultiVectorIterator);

/********************************************************
* Constructors                                          *
********************************************************/
MultiVectorIterator::MultiVectorIterator()
{
    typeID = MultiVectorIteratorTypeID;
    iterator = NULL;
    d_elements = boost::shared_ptr<std::vector<MeshElement> >();
    d_pos = 0;
}
MultiVectorIterator::MultiVectorIterator( boost::shared_ptr<std::vector<MeshElement> > elements, size_t pos )
{
    typeID = MultiVectorIteratorTypeID;
    iterator = NULL;
    d_elements = elements;
    d_pos = pos;
}
MultiVectorIterator::MultiVectorIterator(const MultiVectorIterator& rhs)
{
    typeID = MultiVectorIteratorTypeID;
    iterator = NULL;
    d_elements = rhs.d_elements;
    d_pos = rhs.d_pos;
}
MultiVectorIterator& MultiVectorIterator::operator=(const MultiVectorIterator& rhs)
{
    if (this == &rhs) // protect against invalid self-assignment
        return *this;
    this->typeID = MultiVectorIteratorTypeID;
    this->iterator = NULL;
    this->d_elements = rhs.d_elements;
    this->d_pos = rhs.d_pos;
    return *this;
}


/********************************************************
* Function to clone the iterator                        *
********************************************************/
MeshIterator* MultiVectorIterator::clone() const
{
    return new MultiVectorIterator(*this);
}


/********************************************************
* De-constructor                                        *
********************************************************/
MultiVectorIterator::~MultiVectorIterator()
{
}


/********************************************************
* Return an iterator to the beginning or end            *
********************************************************/
MeshIterator MultiVectorIterator::begin() const
{
    return MultiVectorIterator( d_elements, 0 );
}
MeshIterator MultiVectorIterator::end() const
{
    return MultiVectorIterator( d_elements, d_elements->size() );
}


/********************************************************
* Return the number of elements in the iterator         *
********************************************************/
size_t MultiVectorIterator::size() const
{
    return d_elements->size();
}


/********************************************************
* Increment/Decrement the iterator                      *
********************************************************/
MeshIterator& MultiVectorIterator::operator++()
{
    // Prefix increment (increment and return this)
    d_pos++;
    return *this;
}
MeshIterator MultiVectorIterator::operator++(int)
{
    // Postfix increment (increment and return temporary object)
    MultiVectorIterator tmp(*this);     // Create a temporary variable
    this->operator++();             // apply operator
    return tmp;                     // return temporary result
}
MeshIterator& MultiVectorIterator::operator--()
{
    // Prefix decrement (increment and return this)
    d_pos--;
    return *this;
}
MeshIterator MultiVectorIterator::operator--(int)
{
    // Postfix decrement (increment and return temporary object)
    MultiVectorIterator tmp(*this);      // Create a temporary variable
    --(*this);                      // apply operator
    return tmp;                     // return temporary result
}


/********************************************************
* Compare two iterators                                 *
********************************************************/
bool MultiVectorIterator::operator==(const MeshIterator& rhs)
{
    MultiVectorIterator* rhs2 = NULL;
    MultiVectorIterator* tmp = (MultiVectorIterator*) &rhs;     // Convert rhs to a MultiVectorIterator* so we can access the base class members
    if ( typeid(rhs)==typeid(MultiVectorIterator) ) {
        rhs2 = tmp;     // We can safely cast rhs to a MultiVectorIterator
    } else if ( tmp->typeID==MultiVectorIteratorTypeID ) {
        rhs2 = tmp;     // We can safely cast rhs.iterator to a MultiVectorIterator
    } else if ( ((MultiVectorIterator*)tmp->iterator)->typeID==MultiVectorIteratorTypeID ) {
        rhs2 = (MultiVectorIterator*) tmp->iterator;
    } else {
        AMP_ERROR("Error, comparing a MultiVectorIterator iterator to an unknown iterator");
    }
    if ( rhs2 != NULL ) {
        return d_elements.get()==rhs2->d_elements.get() && d_pos==rhs2->d_pos;
    }
    return false;
}
bool MultiVectorIterator::operator!=(const MeshIterator& rhs)
{
    return !((*this)==rhs);
}


/********************************************************
* Dereference the iterator to get the element           *
********************************************************/
MeshElement& MultiVectorIterator::operator*()
{
    if ( d_pos<0 || d_pos>=d_elements->size() )
        AMP_ERROR("Invalid dereference (iterator is out of range");
    return d_elements->operator[](d_pos);
}
MeshElement* MultiVectorIterator::operator->()
{
    if ( d_pos<0 || d_pos>=d_elements->size() )
        AMP_ERROR("Invalid dereference (iterator is out of range");
    return &(d_elements->operator[](d_pos));
}


} // Mesh namespace
} // AMP namespace

