#include "ampmesh/MeshIterator.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Mesh {


// Create a unique id for this class
static unsigned int MeshIteratorTypeID = TYPE_HASH( MeshIterator );

// unused global variable to prevent compiler warning
static MeshElement nullElement;


/********************************************************
* Constructors                                          *
********************************************************/
MeshIterator::MeshIterator()
{
    typeID   = MeshIteratorTypeID;
    iterator = nullptr;
}
MeshIterator::MeshIterator( const MeshIterator &rhs )
{
    typeID   = MeshIteratorTypeID;
    iterator = nullptr;
    if ( rhs.iterator == nullptr && rhs.typeID == MeshIteratorTypeID ) {
        iterator = nullptr;
    } else if ( rhs.typeID != MeshIteratorTypeID ) {
        iterator = rhs.clone();
    } else {
        iterator = rhs.iterator->clone();
    }
}
MeshIterator &MeshIterator::operator=( const MeshIterator &rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    if ( iterator != nullptr ) {
        // Delete the existing element
        delete iterator;
        iterator = nullptr;
    }
    typeID = MeshIteratorTypeID;
    if ( rhs.iterator == nullptr && rhs.typeID == MeshIteratorTypeID ) {
        iterator = nullptr;
    } else if ( rhs.typeID != MeshIteratorTypeID ) {
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
    if ( iterator != nullptr )
        delete iterator;
    iterator = nullptr;
}


/********************************************************
* Function to clone the iterator                        *
********************************************************/
MeshIterator *MeshIterator::clone() const
{
    if ( iterator == nullptr )
        return new MeshIterator();
    else
        AMP_ERROR( "clone must instantiated by the derived class" );
    return nullptr;
}


/********************************************************
* Functions to return the begin or end iterator         *
********************************************************/
MeshIterator MeshIterator::begin() const
{
    if ( iterator == nullptr )
        return MeshIterator();
    return iterator->begin();
}
MeshIterator MeshIterator::end() const
{
    if ( iterator == nullptr )
        return MeshIterator();
    return iterator->end();
}


/********************************************************
* Functions for incrementing/decrementing               *
********************************************************/
MeshIterator &MeshIterator::operator++() { return iterator->operator++(); }
MeshIterator MeshIterator::operator++( int i ) { return iterator->operator++( i ); }
MeshIterator &MeshIterator::operator--() { return iterator->operator--(); }
MeshIterator MeshIterator::operator--( int i ) { return iterator->operator--( i ); }


/********************************************************
* Functions for incrementing/decrementing               *
********************************************************/
bool MeshIterator::operator==( const MeshIterator &rhs ) const
{
    if ( this->size() == 0 && rhs.size() == 0 )
        return true;
    if ( this->size() != rhs.size() || this->position() != rhs.position() )
        return false;
    if ( iterator == nullptr )
        return rhs.iterator == nullptr;
    return iterator->operator==( rhs );
}
bool MeshIterator::operator!=( const MeshIterator &rhs ) const
{
    if ( this->size() == 0 && rhs.size() == 0 )
        return false;
    if ( this->size() != rhs.size() || this->position() != rhs.position() )
        return true;
    if ( iterator == nullptr )
        return rhs.iterator != nullptr;
    return iterator->operator!=( rhs );
}


/********************************************************
* Functions for dereferencing the iterator              *
********************************************************/
MeshElement &MeshIterator::operator*() { return iterator->operator*(); }
MeshElement *MeshIterator::operator->() { return iterator->operator->(); }
const MeshElement &MeshIterator::operator*() const { return iterator->operator*(); }
const MeshElement *MeshIterator::operator->() const { return iterator->operator->(); }
MeshElement &MeshIterator::operator[]( int i )
{
    if ( iterator != nullptr )
        return iterator->operator[]( i );
    AMP_ERROR( "Dereferencing iterator with offset is not supported by default" );
    return this->operator*(); // This line never executes and would return the wrong object
}


/********************************************************
* Function to get the size and position of the iterator *
********************************************************/
size_t MeshIterator::size() const
{
    if ( iterator == nullptr )
        return 0;
    return iterator->size();
}
size_t MeshIterator::position() const
{
    if ( iterator == nullptr )
        return 0;
    return iterator->position();
}


/********************************************************
*  arithmetic operators                                 *
********************************************************/
MeshIterator MeshIterator::operator+( int n ) const
{
    if ( iterator != nullptr )
        return iterator->operator+( n );
    MeshIterator tmp( *this ); // Create a temporary iterator
    tmp.operator+=( n );       // Increment temporary iterator
    return tmp;                // return temporary iterator
}
MeshIterator MeshIterator::operator+( const MeshIterator &it ) const
{
    if ( iterator != nullptr )
        return iterator->operator+( it );
    return this->operator+( (int) it.position() );
}
MeshIterator MeshIterator::operator-( int n ) const
{
    if ( iterator != nullptr )
        return iterator->operator-( n );
    MeshIterator tmp( *this ); // Create a temporary iterator
    return this->operator+( -n );
}
MeshIterator MeshIterator::operator-( const MeshIterator &it ) const
{
    if ( iterator != nullptr )
        return iterator->operator+( it );
    return this->operator+( -static_cast<int>( it.position() ) );
}
MeshIterator &MeshIterator::operator+=( int n )
{
    if ( iterator != nullptr )
        return iterator->operator+=( n );
    if ( n >= 0 ) {
        for ( int i = 0; i < n; i++ ) {
            this->operator++();
        } // increment iterator
    } else {
        for ( int i = 0; i < -n; i++ ) {
            this->operator--();
        } // decrement iterator
    }
    return *this;
}
MeshIterator &MeshIterator::operator+=( const MeshIterator &it )
{
    if ( iterator != nullptr )
        return iterator->operator+=( (int) it.position() );
    return this->operator+=( (int) it.position() );
}
MeshIterator &MeshIterator::operator-=( int n )
{
    if ( iterator != nullptr )
        return iterator->operator-=( n );
    return this->operator+=( -n );
}
MeshIterator &MeshIterator::operator-=( const MeshIterator &it )
{
    if ( iterator != nullptr )
        return iterator->operator-=( (int) it.position() );
    return this->operator+=( -static_cast<int>( it.position() ) );
}


/********************************************************
*  Comparison operators                                 *
********************************************************/
bool MeshIterator::operator<( const MeshIterator &rhs )
{
    return this->position() < rhs.position();
}
bool MeshIterator::operator<=( const MeshIterator &rhs )
{
    return this->position() <= rhs.position();
}
bool MeshIterator::operator>( const MeshIterator &rhs )
{
    return this->position() > rhs.position();
}
bool MeshIterator::operator>=( const MeshIterator &rhs )
{
    return this->position() >= rhs.position();
}


} // Mesh namespace
} // AMP namespace
