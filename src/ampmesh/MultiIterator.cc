#include "ampmesh/MultiIterator.h"
#include "ampmesh/MeshElement.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Mesh {


// Create a unique id for this class
static unsigned int MultiIteratorTypeID = TYPE_HASH(MultiIterator);

// unused global variable to prevent compiler warning
static MeshElement nullElement;


/********************************************************
* Constructors                                          *
********************************************************/
MultiIterator::MultiIterator()
{
    typeID = MultiIteratorTypeID;
    iterator = NULL;
    d_iterators = std::vector<boost::shared_ptr<MeshIterator> >(0);
    d_iteratorSize = std::vector<size_t>(0);
    d_localPos = 0;
    d_globalPos = 0;
    d_iteratorNum = 0;
}
MultiIterator::MultiIterator( std::vector<boost::shared_ptr<MeshIterator> > iterators, size_t global_pos )
{
    typeID = MultiIteratorTypeID;
    iterator = NULL;
    d_iterators = iterators;
    d_iteratorSize = std::vector<size_t>(d_iterators.size(),0);
    d_globalSize = 0;
    for (size_t i=0; i<d_iterators.size(); i++) {
        d_iteratorSize[i] = d_iterators[i]->size();
        d_globalSize += d_iteratorSize[i];
    }
    d_globalPos = global_pos;
    cur_iterator = MeshIterator();
    // Set the local position and the local iterator    
    if ( d_globalPos > d_globalSize ) {
        // The position is more than one past the last element
        AMP_ERROR("Cannot create a MultiIterator with a current index that is more than 1 past the last point");
    } else if ( d_globalPos == d_globalSize ) {
        // The position is one past the last element
        d_localPos = 0;
        d_iteratorNum = d_iterators.size();
    } else {
        // We are inside the iterator, we need to point to the current element
        d_iteratorNum = 0;
        d_localPos = global_pos;
        while ( d_localPos >= d_iteratorSize[d_iteratorNum] ) {
            d_iteratorNum++;
            d_localPos -= d_iteratorSize[d_iteratorNum];
        }
        cur_iterator = d_iterators[d_iteratorNum]->begin();
        for (size_t i=0; i<d_localPos; i++)
            ++cur_iterator;
    }
}
MultiIterator::MultiIterator(const MultiIterator& rhs)
{
    typeID = MultiIteratorTypeID;
    iterator = NULL;
    d_iterators    = rhs.d_iterators;
    d_iteratorSize = rhs.d_iteratorSize;
    d_globalSize   = rhs.d_globalSize;
    d_localPos     = rhs.d_localPos;
    d_globalPos    = rhs.d_globalPos;
    d_iteratorNum  = rhs.d_iteratorNum;
    cur_iterator   = rhs.cur_iterator;
}
MultiIterator& MultiIterator::operator=(const MultiIterator& rhs)
{
    if (this == &rhs) // protect against invalid self-assignment
        return *this;
    this->typeID = MultiIteratorTypeID;
    this->iterator = NULL;
    this->d_iterators    = rhs.d_iterators;
    this->d_iteratorSize = rhs.d_iteratorSize;
    this->d_globalSize   = rhs.d_globalSize;
    this->d_localPos     = rhs.d_localPos;
    this->d_globalPos    = rhs.d_globalPos;
    this->d_iteratorNum  = rhs.d_iteratorNum;
    this->cur_iterator   = rhs.cur_iterator;
    return *this;
}


/********************************************************
* Function to clone the iterator                        *
********************************************************/
MeshIterator* MultiIterator::clone() const
{
    return new MultiIterator(*this);
}


/********************************************************
* De-constructor                                        *
********************************************************/
MultiIterator::~MultiIterator()
{
}


/********************************************************
* Return an iterator to the beginning or end            *
********************************************************/
MeshIterator MultiIterator::begin() const
{
    return MultiIterator( d_iterators, 0 );
}
MeshIterator MultiIterator::end() const
{
    return MultiIterator( d_iterators, d_globalSize );
}


/********************************************************
* Return the number of elements in the iterator         *
********************************************************/
size_t MultiIterator::size() const
{
    return d_globalSize;
}


/********************************************************
* Increment/Decrement the iterator                      *
********************************************************/
MeshIterator& MultiIterator::operator++()
{
    // Prefix increment (increment and return this)
    if ( d_globalPos == d_globalSize )
        AMP_ERROR("Iterating more than ont past the last element");
    if ( d_globalPos+1 == d_globalSize ) {
        // We have moved to one past the last element
        d_globalPos = d_globalSize;
        d_localPos = 0;
        d_iteratorNum = d_iterators.size();
        cur_iterator = MeshIterator();
    } else if ( d_localPos+1 == d_iteratorSize[d_iteratorNum] ) {
        // We need to change the internal iterator
        d_globalPos++;
        d_localPos = 0;
        d_iteratorNum++;
        cur_iterator = d_iterators[d_iteratorNum]->begin();
    } else {
        // We are within the same iterator
        d_localPos++;
        d_globalPos++;
        ++cur_iterator;     // preincrement for consistency and speed
    }
    return *this;
}
MeshIterator MultiIterator::operator++(int)
{
    // Postfix increment (increment and return temporary object)
    MultiIterator tmp(*this);     // Create a temporary variable
    this->operator++();             // apply operator
    return tmp;                     // return temporary result
}
MeshIterator& MultiIterator::operator--()
{
    // Prefix decrement (increment and return this)
    AMP_ERROR("Decrementing MultiMesh iterators is not implimented yet");
    return *this;
}
MeshIterator MultiIterator::operator--(int)
{
    // Postfix decrement (increment and return temporary object)
    MultiIterator tmp(*this);      // Create a temporary variable
    --(*this);                      // apply operator
    return tmp;                     // return temporary result
}


/********************************************************
* Compare two iterators                                 *
* Two MultiIterators are the same if both the list of   *
* iterators and the current position are the same.      *
********************************************************/
bool MultiIterator::operator==(const MeshIterator& rhs)
{
    MultiIterator* rhs2 = NULL;
    MultiIterator* tmp = (MultiIterator*) &rhs;     // Convert rhs to a MultiIterator* so we can access the base class members
    if ( typeid(rhs)==typeid(MultiIterator) ) {
        rhs2 = tmp;     // We can safely cast rhs to a MultiIterator
    } else if ( tmp->typeID==MultiIteratorTypeID ) {
        rhs2 = tmp;     // We can safely cast rhs.iterator to a MultiIterator
    } else if ( ((MultiIterator*)tmp->iterator)->typeID==MultiIteratorTypeID ) {
        rhs2 = (MultiIterator*) tmp->iterator;
    } else {
        AMP_ERROR("Error, comparing a MultiIterator iterator to an unknown iterator");
    }
    if ( rhs2 != NULL ) {
        bool equal = true;
        equal = equal && d_globalSize==rhs2->d_globalSize;
        equal = equal && d_globalPos==rhs2->d_globalPos;
        equal = equal && d_iterators.size()==rhs2->d_iterators.size();
        if ( equal ) {
            for (size_t i=0; i<d_iterators.size(); i++)
                equal = equal && (*d_iterators[i])==(*rhs2->d_iterators[i]);
        }
        return equal;
    }
    return false;
}
bool MultiIterator::operator!=(const MeshIterator& rhs)
{
    return !((*this)==rhs);
}


/********************************************************
* Dereference the iterator to get the element           *
********************************************************/
MeshElement& MultiIterator::operator*()
{
    if ( d_globalPos >= d_globalSize )
        AMP_ERROR("Invalid dereference (iterator is out of range");
    return cur_iterator.operator*();
}
MeshElement* MultiIterator::operator->()
{
    if ( d_globalPos >= d_globalSize )
        AMP_ERROR("Invalid dereference (iterator is out of range");
    return cur_iterator.operator->();
}


} // Mesh namespace
} // AMP namespace

