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
* Functions that aren't implimented for the base class  *
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

