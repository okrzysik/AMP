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
    d_iterators = std::vector<AMP::shared_ptr<MeshIterator> >(0);
    d_iteratorSize = std::vector<size_t>(0);
    d_localPos = 0;
    d_globalPos = 0;
    d_iteratorNum = 0;
}
MultiIterator::MultiIterator( std::vector<AMP::shared_ptr<MeshIterator> > iterators, size_t global_pos )
{
    typeID = MultiIteratorTypeID;
    iterator = NULL;
    d_iterators.resize(0);
    for (size_t i=0; i<iterators.size(); i++) {
        if ( iterators[i]==NULL )
            continue;
        if ( iterators[i]->size() > 0 )
            d_iterators.push_back(iterators[i]);
    }
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
size_t MultiIterator::position() const
{
    return d_globalPos;
}


/********************************************************
* Increment/Decrement the iterator                      *
********************************************************/
MeshIterator& MultiIterator::operator++()
{
    // Prefix increment (increment and return this)
    if ( d_globalPos == d_globalSize )
        AMP_ERROR("Iterating more than one past the last element");
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
    if ( d_globalPos == 0 )
        AMP_ERROR("Iterating before the first element");
    if ( d_globalPos == d_globalSize ) {
        // We are starting at the end
        d_globalPos = d_globalSize+1;
        d_iteratorNum = d_iterators.size()-1;
        d_localPos = d_iteratorSize[d_iteratorNum]-1;
        cur_iterator = d_iterators[d_iteratorNum]->end();
        --cur_iterator;
    } else if ( d_localPos == 0 ) {
        // We need to change the internal iterator
        d_globalPos--;
        d_iteratorNum--;
        d_localPos = d_iteratorSize[d_iteratorNum]-1;
        cur_iterator = d_iterators[d_iteratorNum]->end();
        --cur_iterator;
    } else {
        // We are within the same iterator
        d_localPos--;
        d_globalPos--;
        --cur_iterator;     // predecrement for consistency and speed
    }
    return *this;
}
MeshIterator MultiIterator::operator--(int)
{
    // Postfix decrement (increment and return temporary object)
    MultiIterator tmp(*this);       // Create a temporary variable
    --(*this);                      // apply operator
    return tmp;                     // return temporary result
}


/********************************************************
* Random access incrementors                            *
********************************************************/
MeshIterator MultiIterator::operator+(int n) const
{
    MultiIterator tmp(*this);           // Create a temporary iterator
    tmp.operator+=(n);                  // Increment temporary iterator
    return tmp;
}
MeshIterator& MultiIterator::operator+=(int n)
{
    if ( n>=0 ) {                       // increment *this
        size_t n2 = static_cast<size_t>(n);
        if ( d_globalPos+n2 > d_globalSize )
            AMP_ERROR("Iterated past end of iterator");
        if ( d_globalPos+n2 == d_globalSize ) {
            // We reached the end of the iterator
            d_globalPos = d_globalSize;
            d_localPos = 0;
            d_iteratorNum = d_iterators.size();
            cur_iterator = MeshIterator();
            return *this;
        }
        // Move to the correct iterator
        if ( d_localPos+n2 >= d_iteratorSize[d_iteratorNum] ) {
            size_t i = d_iteratorSize[d_iteratorNum] - d_localPos;
            n2 -= i;
            d_iteratorNum++;
            while ( d_iteratorSize[d_iteratorNum] < n2 ) {
                n2 -= d_iteratorSize[d_iteratorNum];
                d_iteratorNum++;
            }
            cur_iterator = d_iterators[d_iteratorNum]->begin();
            d_localPos = 0;
        }
        // Increment local iterator
        cur_iterator.operator+=(n2);
        d_localPos += n2;
        d_globalPos += n;
    } else {                            // decrement *this
        size_t n2 = static_cast<size_t>(-n);
        if ( d_globalPos < n2 )
            AMP_ERROR("Iterated past beginning of iterator");
        if ( d_globalPos == n2 ) {
            // We reached the beginning of the iterator
            d_globalPos = 0;
            d_iteratorNum = 0;
            d_localPos = 0;
            cur_iterator = d_iterators[0]->begin();
            return *this;
        }
        // Move to the correct iterator
        if ( n2 > d_localPos ) {
            //size_t i = d_iteratorSize[d_iteratorNum] - d_localPos;
            n2 -= d_localPos;
            d_iteratorNum--;
            while ( d_iteratorSize[d_iteratorNum] < n2 ) {
                n2 -= d_iteratorSize[d_iteratorNum];
                d_iteratorNum--;
            }
            cur_iterator = d_iterators[d_iteratorNum]->end();
            d_localPos = d_iteratorSize[d_iteratorNum];
        }
        // Increment local iterator
        cur_iterator.operator-=(n2);
        d_localPos -= n2;
        d_globalPos += n;
    }
    return *this;
}


/********************************************************
* Compare two iterators                                 *
* Two MultiIterators are the same if both the list of   *
* iterators and the current position are the same.      *
********************************************************/
bool MultiIterator::operator==(const MeshIterator& rhs) const
{
    MultiIterator* rhs2 = NULL;
    MultiIterator* tmp = (MultiIterator*) &rhs;     // Convert rhs to a MultiIterator* so we can access the base class members
    if ( typeid(rhs)==typeid(MultiIterator) ) {
        rhs2 = tmp;     // We can safely cast rhs to a MultiIterator
    } else if ( tmp->typeID==MultiIteratorTypeID ) {
        rhs2 = tmp;     // We can safely cast rhs.iterator to a MultiIterator
    } else if ( ((MultiIterator*)tmp->iterator)->typeID==MultiIteratorTypeID ) {
        rhs2 = (MultiIterator*) tmp->iterator;
    }
    // Perform direct comparisions if we are dealing with two MultiIterator
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
    /* We are comparing a MultiIterator to an arbitrary iterator
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
    MeshIterator iterator1 = this->begin();
    MeshIterator iterator2 = rhs.begin();
    bool elements_match = true;
    for (size_t i=0; i<this->size(); i++) {
        if ( iterator1->globalID() != iterator2->globalID() )
            elements_match = false;
        ++iterator1;
        ++iterator2;
    }
    return elements_match;
}
bool MultiIterator::operator!=(const MeshIterator& rhs) const
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

