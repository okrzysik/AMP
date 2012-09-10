#include "ampmesh/MeshIterator.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Mesh {


// Create a unique id for this class
static unsigned int MeshIteratorTypeID = TYPE_HASH(MeshIterator);

// unused global variable to prevent compiler warning
static MeshElement nullElement;



/********************************************************
* Constructors                                          *
********************************************************/
MeshIterator::MeshIterator()
{
    typeID = MeshIteratorTypeID;
    iterator = NULL;
}
MeshIterator::MeshIterator(const MeshIterator& rhs)
{
    typeID = MeshIteratorTypeID;
    iterator = NULL;
    if ( rhs.iterator==NULL && rhs.typeID==MeshIteratorTypeID ) {
        iterator = NULL;
    } else if ( rhs.typeID!=MeshIteratorTypeID ) {
        iterator = rhs.clone();
    } else {
        iterator = rhs.iterator->clone();
    }
}
MeshIterator& MeshIterator::operator=(const MeshIterator& rhs)
{
    if (this == &rhs) // protect against invalid self-assignment
        return *this;
    if ( iterator != NULL ) {
        // Delete the existing element
        delete iterator;
        iterator = NULL;
    }
    typeID = MeshIteratorTypeID;
    if ( rhs.iterator==NULL && rhs.typeID==MeshIteratorTypeID ) {
        iterator = NULL;
    } else if ( rhs.typeID!=MeshIteratorTypeID ) {
        iterator = rhs.clone();
    } else {
        iterator = rhs.iterator->clone();
    }
    return *this;
}


/********************************************************
* De-constructor                                        *
********************************************************/
MeshIterator::~MeshIterator()
{
    if ( iterator != NULL )
        delete iterator;
    iterator = NULL;
}


/********************************************************
* Function to clone the iterator                        *
********************************************************/
MeshIterator* MeshIterator::clone() const
{
    if ( iterator==NULL )
        return new MeshIterator();
    else
        AMP_ERROR("clone must instantiated by the derived class");
    return NULL;
}


/********************************************************
* Functions to return the begin or end iterator         *
********************************************************/
MeshIterator MeshIterator::begin() const
{
    if ( iterator==NULL )
        return MeshIterator();
    return iterator->begin();
}
MeshIterator MeshIterator::end() const
{
    if ( iterator==NULL )
        return MeshIterator();
    return iterator->end();
}


/********************************************************
* Functions for incrementing/decrementing               *
********************************************************/
MeshIterator& MeshIterator::operator++()
{
    return iterator->operator++();
}
MeshIterator MeshIterator::operator++(int i)
{
    return iterator->operator++(i);
}
MeshIterator& MeshIterator::operator--()
{
    return iterator->operator--();
}
MeshIterator MeshIterator::operator--(int i)
{
    return iterator->operator--(i);
}


/********************************************************
* Functions for incrementing/decrementing               *
********************************************************/
bool MeshIterator::operator==(const MeshIterator& rhs) const
{
    if ( iterator==NULL )
        return rhs.iterator==NULL;
    return iterator->operator==(rhs);
}
bool MeshIterator::operator!=(const MeshIterator& rhs) const
{
    if ( iterator==NULL )
        return rhs.iterator!=NULL;
    return iterator->operator!=(rhs);
}


/********************************************************
* Functions for dereferencing the iterator              *
********************************************************/
MeshElement& MeshIterator::operator*()
{
    return iterator->operator*();
}
MeshElement* MeshIterator::operator->()
{
    return iterator->operator->();
}
MeshElement& MeshIterator::operator[](int i)
{
    if ( iterator!=NULL )
        return iterator->operator[](i);
    AMP_ERROR("Dereferencing iterator with offset is not supported by default");
    return this->operator*();   // This line never executes and would return the wrong object
}


/********************************************************
* Function to get the size and position of the iterator *
********************************************************/
size_t MeshIterator::size() const
{
    if ( iterator==NULL )
        return 0;
    return iterator->size();
}
size_t MeshIterator::position() const
{
    if ( iterator==NULL )
        return 0;
    return iterator->position();
}


/********************************************************
*  arithmetic operators                                 *
********************************************************/
MeshIterator MeshIterator::operator+(int n) const
{
    MeshIterator tmp(*this);            // Create a temporary iterator
    for (int i=0; i<n; i++) { ++tmp; }  // increment temporary iterator
    return tmp;                         // return temporary iterator
}
MeshIterator MeshIterator::operator+(const MeshIterator& it) const
{
    return this->operator+((int)it.position());
}
MeshIterator MeshIterator::operator-(int n) const
{
    MeshIterator tmp(*this);            // Create a temporary iterator
    for (int i=0; i<n; i++) { --tmp; }  // decrement temporary iterator
    return tmp;                         // return temporary iterator
}
MeshIterator MeshIterator::operator-(const MeshIterator& it) const
{
    return this->operator-((int)it.position());
}
MeshIterator& MeshIterator::operator+=(int n)
{
    if ( iterator!=NULL ) 
        return iterator->operator+=(n);
    for (int i=0; i<n; i++)
        this->operator++();
    return *this;
}
MeshIterator& MeshIterator::operator+=(const MeshIterator& it)
{
    if ( iterator!=NULL )
        return iterator->operator+=((int)it.position());
    return this->operator+=((int)it.position());
}
MeshIterator& MeshIterator::operator-=(int n)
{
    if ( iterator!=NULL )
        return iterator->operator-=(n);
    for (int i=0; i<n; i++)
        this->operator--();
    return *this;
}
MeshIterator& MeshIterator::operator-=(const MeshIterator& it)
{
    if ( iterator!=NULL )
        return iterator->operator-=((int)it.position());
    return this->operator-=((int)it.position());
}


/********************************************************
*  Comparison operators                                 *
********************************************************/
bool MeshIterator::operator<(const MeshIterator& rhs)
{
    return this->position() < rhs.position();
}
bool MeshIterator::operator<=(const MeshIterator& rhs)
{
    return this->position() <= rhs.position();
}
bool MeshIterator::operator>(const MeshIterator& rhs)
{
    return this->position() > rhs.position();
}
bool MeshIterator::operator>=(const MeshIterator& rhs)
{
    return this->position() >= rhs.position();
}


} // Mesh namespace
} // AMP namespace

