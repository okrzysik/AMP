#include "AMP/ampmesh/MeshIterator.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Mesh {


// unused global variable to prevent compiler warning
static MeshElement nullElement;


/********************************************************
 * Set the base class type id                            *
 ********************************************************/
const uint32_t MeshIterator::MeshIteratorTypeID = AMP::Utilities::hash_char( "MeshIterator" );


/********************************************************
 * Constructors                                          *
 ********************************************************/
MeshIterator::MeshIterator()
    : d_iterator( nullptr ),
      d_typeID( MeshIteratorTypeID ),
      d_iteratorType( Type::RandomAccess ),
      d_size( 0 ),
      d_pos( 0 ),
      d_element( nullptr )
{
}
MeshIterator::MeshIterator( MeshIterator &&rhs )
    : d_iterator( nullptr ),
      d_typeID( MeshIteratorTypeID ),
      d_iteratorType( rhs.d_iteratorType ),
      d_size( 0 ),
      d_pos( 0 ),
      d_element( nullptr )
{
    if ( rhs.d_iterator == nullptr && rhs.d_typeID == MeshIteratorTypeID ) {
        d_iterator = nullptr;
    } else if ( rhs.d_typeID != MeshIteratorTypeID ) {
        d_iterator = rhs.clone();
    } else {
        d_iterator     = rhs.d_iterator;
        rhs.d_iterator = nullptr;
    }
}
MeshIterator::MeshIterator( const MeshIterator &rhs )
    : d_iterator( nullptr ),
      d_typeID( MeshIteratorTypeID ),
      d_iteratorType( rhs.d_iteratorType ),
      d_size( 0 ),
      d_pos( 0 ),
      d_element( nullptr )
{
    if ( rhs.d_iterator == nullptr && rhs.d_typeID == MeshIteratorTypeID ) {
        d_iterator = nullptr;
    } else if ( rhs.d_typeID != MeshIteratorTypeID ) {
        d_iterator = rhs.clone();
    } else {
        d_iterator = rhs.d_iterator->clone();
    }
}
MeshIterator &MeshIterator::operator=( MeshIterator &&rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    if ( d_iterator != nullptr ) {
        // Delete the existing element
        delete d_iterator;
        d_iterator = nullptr;
    }
    d_typeID       = MeshIteratorTypeID;
    d_iteratorType = rhs.d_iteratorType;
    d_size         = 0;
    d_pos          = 0;
    d_element      = nullptr;
    if ( rhs.d_iterator == nullptr && rhs.d_typeID == MeshIteratorTypeID ) {
        d_iterator = nullptr;
    } else if ( rhs.d_typeID != MeshIteratorTypeID ) {
        d_iterator = rhs.clone();
    } else {
        d_iterator     = rhs.d_iterator;
        rhs.d_iterator = nullptr;
    }
    return *this;
}
MeshIterator &MeshIterator::operator=( const MeshIterator &rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    if ( d_iterator != nullptr ) {
        // Delete the existing element
        delete d_iterator;
        d_iterator = nullptr;
    }
    d_typeID       = MeshIteratorTypeID;
    d_iteratorType = rhs.d_iteratorType;
    d_size         = 0;
    d_pos          = 0;
    d_element      = nullptr;
    if ( rhs.d_iterator == nullptr && rhs.d_typeID == MeshIteratorTypeID ) {
        d_iterator = nullptr;
    } else if ( rhs.d_typeID != MeshIteratorTypeID ) {
        d_iterator = rhs.clone();
    } else {
        d_iterator = rhs.d_iterator->clone();
    }
    return *this;
}


/********************************************************
 * De-constructor                                        *
 ********************************************************/
MeshIterator::~MeshIterator()
{
    if ( d_iterator != nullptr )
        delete d_iterator;
    d_iterator = nullptr;
}


/********************************************************
 * Function to clone the d_iterator                      *
 ********************************************************/
MeshIterator *MeshIterator::clone() const
{
    if ( d_iterator == nullptr )
        return new MeshIterator();
    else
        AMP_ERROR( "clone must be instantiated by the derived class" );
    return nullptr;
}


/********************************************************
 * Return the iterator type                              *
 ********************************************************/
MeshIterator::Type MeshIterator::type() const
{
    if ( d_iterator == nullptr )
        return d_iteratorType;
    return d_iterator->d_iteratorType;
}


/********************************************************
 * Functions to return the begin or end d_iterator       *
 ********************************************************/
MeshIterator MeshIterator::begin() const
{
    if ( d_iterator == nullptr )
        return MeshIterator();
    return d_iterator->begin();
}
MeshIterator MeshIterator::end() const
{
    if ( d_iterator == nullptr )
        return MeshIterator();
    return d_iterator->end();
}


/********************************************************
 * Functions for incrementing/decrementing               *
 ********************************************************/
MeshIterator &MeshIterator::operator++() { return d_iterator->operator++(); }
MeshIterator &MeshIterator::operator--() { return d_iterator->operator--(); }
MeshIterator MeshIterator::operator++( int i ) { return d_iterator->operator++( i ); }
MeshIterator MeshIterator::operator--( int i ) { return d_iterator->operator--( i ); }


/********************************************************
 * Functions for incrementing/decrementing               *
 ********************************************************/
bool MeshIterator::operator==( const MeshIterator &rhs ) const
{
    if ( this->size() == 0 && rhs.size() == 0 )
        return true;
    if ( this->size() != rhs.size() || this->position() != rhs.position() )
        return false;
    if ( d_iterator == nullptr )
        return rhs.d_iterator == nullptr;
    return d_iterator->operator==( rhs );
}
bool MeshIterator::operator!=( const MeshIterator &rhs ) const
{
    if ( this->size() == 0 && rhs.size() == 0 )
        return false;
    if ( this->size() != rhs.size() || this->position() != rhs.position() )
        return true;
    if ( d_iterator == nullptr )
        return rhs.d_iterator != nullptr;
    return d_iterator->operator!=( rhs );
}


/********************************************************
 * Functions for dereferencing the d_iterator            *
 ********************************************************/
MeshElement &MeshIterator::operator[]( int i )
{
    if ( d_iterator != nullptr )
        return d_iterator->operator[]( i );
    AMP_ERROR( "Dereferencing d_iterator with offset is not supported by default" );
    return this->operator*(); // This line never executes and would return the wrong object
}


/********************************************************
 *  arithmetic operators                                 *
 ********************************************************/
MeshIterator MeshIterator::operator+( int n ) const
{
    if ( d_iterator != nullptr )
        return d_iterator->operator+( n );
    MeshIterator tmp( *this ); // Create a temporary d_iterator
    tmp.operator+=( n );       // Increment temporary d_iterator
    return tmp;                // return temporary d_iterator
}
MeshIterator MeshIterator::operator+( const MeshIterator &it ) const
{
    if ( d_iterator != nullptr )
        return d_iterator->operator+( it );
    return this->operator+( (int) it.position() );
}
MeshIterator MeshIterator::operator-( int n ) const
{
    if ( d_iterator != nullptr )
        return d_iterator->operator-( n );
    MeshIterator tmp( *this ); // Create a temporary d_iterator
    return this->operator+( -n );
}
MeshIterator MeshIterator::operator-( const MeshIterator &it ) const
{
    if ( d_iterator != nullptr )
        return d_iterator->operator+( it );
    return this->operator+( -static_cast<int>( it.position() ) );
}
MeshIterator &MeshIterator::operator+=( int n )
{
    if ( d_iterator != nullptr )
        return d_iterator->operator+=( n );
    if ( n >= 0 ) {
        for ( int i = 0; i < n; i++ ) {
            this->operator++();
        } // increment d_iterator
    } else {
        for ( int i = 0; i < -n; i++ ) {
            this->operator--();
        } // decrement d_iterator
    }
    return *this;
}
MeshIterator &MeshIterator::operator+=( const MeshIterator &it )
{
    if ( d_iterator != nullptr )
        return d_iterator->operator+=( (int) it.position() );
    return this->operator+=( (int) it.position() );
}
MeshIterator &MeshIterator::operator-=( int n )
{
    if ( d_iterator != nullptr )
        return d_iterator->operator-=( n );
    return this->operator+=( -n );
}
MeshIterator &MeshIterator::operator-=( const MeshIterator &it )
{
    if ( d_iterator != nullptr )
        return d_iterator->operator-=( (int) it.position() );
    return this->operator+=( -static_cast<int>( it.position() ) );
}


} // namespace Mesh
} // namespace AMP
