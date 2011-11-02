#include "ampmesh/MeshIterator.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Mesh {


// Create a unique id for this class
static int MeshIteratorTypeID = 1;

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
    if ( rhs.iterator==NULL ) {
        iterator = NULL;
    } else {
        iterator = rhs.clone();
    }
}
MeshIterator& MeshIterator::operator=(const MeshIterator& rhs)
{
    typeID = MeshIteratorTypeID;
    if (this == &rhs) // protect against invalid self-assignment
        return *this;
    if ( rhs.iterator==NULL ) {
        iterator = NULL;
    } else {
        iterator = rhs.clone();
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
* Functions that aren't implimented for the base class  *
********************************************************/
MeshIterator& MeshIterator::operator++()
{
    (*iterator)++;
    return *this;
}
MeshIterator MeshIterator::operator++(int)
{
    return (*iterator)++;
}
MeshIterator& MeshIterator::operator--()
{
    (*iterator)--;
    return *this;
}
MeshIterator MeshIterator::operator--(int)
{
    return (*iterator)++;
}
bool MeshIterator::operator==(const MeshIterator& rhs)
{
    return iterator->operator==(rhs);
}
bool MeshIterator::operator!=(const MeshIterator& rhs)
{
    return iterator->operator!=(rhs);
}
MeshElement& MeshIterator::operator*()
{
    return iterator->operator*();
}
MeshElement* MeshIterator::operator->()
{
    return iterator->operator->();
}
MeshIterator MeshIterator::begin()
{
    return iterator->begin();
}
MeshIterator MeshIterator::end()
{
    return iterator->end();
}



} // Mesh namespace
} // AMP namespace

