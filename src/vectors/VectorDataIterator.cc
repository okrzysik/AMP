#include "vectors/VectorDataIterator.h"
#include "vectors/Vector.h"

namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* Contructors                                                   *
****************************************************************/
VectorDataIterator::VectorDataIterator()
{
    d_Vec          = nullptr;
    d_Block        = nullptr;
    d_CurBlock     = 0;
    d_CurOffset    = 0;
    d_position     = 0;
    d_size         = 0;
    d_CurBlockSize = 0;
}
VectorDataIterator::VectorDataIterator( const VectorDataIterator &rhs )
{
    d_Vec          = rhs.d_Vec;
    d_Block        = rhs.d_Block;
    d_CurBlock     = rhs.d_CurBlock;
    d_CurOffset    = rhs.d_CurOffset;
    d_position     = rhs.d_position;
    d_size         = rhs.d_size;
    d_CurBlockSize = rhs.d_CurBlockSize;
}
VectorDataIterator::VectorDataIterator( Vector *p, size_t position )
{
    d_Vec          = p;
    d_Block        = nullptr;
    d_CurBlock     = 0;
    d_CurOffset    = 0;
    d_position     = 0;
    d_size         = d_Vec->getLocalSize();
    d_CurBlockSize = 0;
    while ( d_position < position ) {
        if ( position - d_position >= d_Vec->sizeOfDataBlock( d_CurBlock ) ) {
            // We need to advance to the next data block
            d_position += d_Vec->sizeOfDataBlock( d_CurBlock );
            d_CurBlock++;
        } else {
            // We need to advance to the correct position within the current data block
            d_CurOffset = position - d_position;
            d_position += d_CurOffset;
        }
    }
    if ( d_CurBlock < d_Vec->numberOfDataBlocks() ) {
        d_Block        = d_Vec->getRawDataBlock<double>( d_CurBlock );
        d_CurBlockSize = d_Vec->sizeOfDataBlock( d_CurBlock );
    }
}
ConstVectorDataIterator::ConstVectorDataIterator()
{
    d_Vec          = nullptr;
    d_Block        = nullptr;
    d_CurBlock     = 0;
    d_CurOffset    = 0;
    d_position     = 0;
    d_size         = 0;
    d_CurBlockSize = 0;
}
ConstVectorDataIterator::ConstVectorDataIterator( const VectorDataIterator &rhs )
{
    d_Vec          = rhs.d_Vec;
    d_Block        = rhs.d_Block;
    d_CurBlock     = rhs.d_CurBlock;
    d_CurOffset    = rhs.d_CurOffset;
    d_position     = rhs.d_position;
    d_size         = rhs.d_size;
    d_CurBlockSize = rhs.d_CurBlockSize;
}
ConstVectorDataIterator::ConstVectorDataIterator( const ConstVectorDataIterator &rhs )
{
    d_Vec          = rhs.d_Vec;
    d_Block        = rhs.d_Block;
    d_CurBlock     = rhs.d_CurBlock;
    d_CurOffset    = rhs.d_CurOffset;
    d_position     = rhs.d_position;
    d_size         = rhs.d_size;
    d_CurBlockSize = rhs.d_CurBlockSize;
}
ConstVectorDataIterator::ConstVectorDataIterator( const Vector *p, size_t position )
{
    d_Vec          = p;
    d_Block        = nullptr;
    d_CurBlock     = 0;
    d_CurOffset    = 0;
    d_position     = 0;
    d_size         = d_Vec->getLocalSize();
    d_CurBlockSize = 0;
    while ( d_position < position ) {
        if ( position - d_position >= d_Vec->sizeOfDataBlock( d_CurBlock ) ) {
            // We need to advance to the next data block
            d_position += d_Vec->sizeOfDataBlock( d_CurBlock );
            d_CurBlock++;
        } else {
            // We need to advance to the correct position within the current data block
            d_CurOffset = position - d_position;
            d_position += d_CurOffset;
        }
    }
    if ( d_CurBlock < d_Vec->numberOfDataBlocks() ) {
        d_Block        = d_Vec->getRawDataBlock<double>( d_CurBlock );
        d_CurBlockSize = d_Vec->sizeOfDataBlock( d_CurBlock );
    }
}


/****************************************************************
* Function to advance/receed the iterator                       *
****************************************************************/
inline void VectorDataIterator::advance( size_t i )
{
    AMP_INSIST( d_position + i <= d_size, "Attempted to iterate past the end of the iterator" );
    size_t togo = i;
    while ( togo > 0 ) {
        if ( togo >= d_CurBlockSize - d_CurOffset ) {
            // We need to advance to the next data block
            d_position += d_CurBlockSize - d_CurOffset;
            togo -= d_CurBlockSize - d_CurOffset;
            d_CurOffset = 0;
            d_CurBlock++;
            if ( d_position < d_size ) {
                d_Block        = d_Vec->getRawDataBlock<double>( d_CurBlock );
                d_CurBlockSize = d_Vec->sizeOfDataBlock( d_CurBlock );
            }
        } else {
            // We need to advance to the correct position within the current data block
            d_position += togo;
            d_CurOffset += togo;
            togo = 0;
        }
    }
}
inline void ConstVectorDataIterator::advance( size_t i )
{
    AMP_INSIST( d_position + i <= d_size, "Attempted to iterate past the end of the iterator" );
    size_t togo = i;
    while ( togo > 0 ) {
        if ( togo >= d_CurBlockSize - d_CurOffset ) {
            // We need to advance to the next data block
            AMP_ASSERT( d_CurBlockSize > 0 );
            d_position += d_CurBlockSize - d_CurOffset;
            togo -= d_CurBlockSize - d_CurOffset;
            d_CurOffset = 0;
            d_CurBlock++;
            if ( d_position < d_size ) {
                d_Block        = d_Vec->getRawDataBlock<double>( d_CurBlock );
                d_CurBlockSize = d_Vec->sizeOfDataBlock( d_CurBlock );
            }
        } else {
            // We need to advance to the correct position within the current data block
            d_position += togo;
            d_CurOffset += togo;
            togo = 0;
        }
    }
}
inline void VectorDataIterator::recede( size_t i )
{
    AMP_INSIST( d_position >= i, "Attempted to iterate past the beginning of the iterator" );
    size_t togo = i;
    while ( togo > 0 ) {
        if ( togo > d_CurOffset ) {
            // We need to advance to the next data block
            d_position -= d_CurOffset;
            togo -= d_CurOffset;
            d_CurBlock--;
            d_Block        = d_Vec->getRawDataBlock<double>( d_CurBlock );
            d_CurBlockSize = d_Vec->sizeOfDataBlock( d_CurBlock );
            d_CurOffset    = d_CurBlockSize;
        } else {
            // We need to advance to the correct position within the current data block
            d_position -= togo;
            d_CurOffset -= togo;
            togo = 0;
        }
    }
}
inline void ConstVectorDataIterator::recede( size_t i )
{
    AMP_INSIST( d_position >= i, "Attempted to iterate past the beginning of the iterator" );
    size_t togo = i;
    while ( togo > 0 ) {
        if ( togo > d_CurOffset ) {
            // We need to advance to the next data block
            d_position -= d_CurOffset;
            togo -= d_CurOffset;
            d_CurBlock--;
            d_Block        = d_Vec->getRawDataBlock<double>( d_CurBlock );
            d_CurBlockSize = d_Vec->sizeOfDataBlock( d_CurBlock );
            d_CurOffset    = d_CurBlockSize;
        } else {
            // We need to advance to the correct position within the current data block
            d_position -= togo;
            d_CurOffset -= togo;
            togo = 0;
        }
    }
}


/****************************************************************
* Increment/Decrement Operators                                 *
****************************************************************/
VectorDataIterator &VectorDataIterator::operator++()
{
    // Prefix increment (increment and return this)
    advance( 1 );
    return *this;
}
VectorDataIterator VectorDataIterator::operator++( int )
{
    // Postfix increment (increment and return temporary object)
    VectorDataIterator tmp( *this ); // Create a temporary variable
    advance( 1 );                    // apply operator
    return tmp;                      // return temporary result
}
ConstVectorDataIterator &ConstVectorDataIterator::operator++()
{
    // Prefix increment (increment and return this)
    advance( 1 );
    return *this;
}
ConstVectorDataIterator ConstVectorDataIterator::operator++( int )
{
    // Postfix increment (increment and return temporary object)
    ConstVectorDataIterator tmp( *this ); // Create a temporary variable
    advance( 1 );                         // apply operator
    return tmp;                           // return temporary result
}
VectorDataIterator &VectorDataIterator::operator--()
{
    // Prefix decrement (decrement and return this)
    recede( 1 );
    return *this;
}
VectorDataIterator VectorDataIterator::operator--( int )
{
    // Postfix decrement (decrement and return temporary object)
    VectorDataIterator tmp( *this ); // Create a temporary variable
    recede( 1 );                     // apply operator
    return tmp;                      // return temporary result
}
ConstVectorDataIterator &ConstVectorDataIterator::operator--()
{
    // Prefix decrement (decrement and return this)
    recede( 1 );
    return *this;
}
ConstVectorDataIterator ConstVectorDataIterator::operator--( int )
{
    // Postfix decrement (decrement and return temporary object)
    ConstVectorDataIterator tmp( *this ); // Create a temporary variable
    recede( 1 );                          // apply operator
    return tmp;                           // return temporary result
}
VectorDataIterator &VectorDataIterator::operator+=( int offset )
{
    if ( offset > 0 )
        advance( offset );
    if ( offset < 0 )
        recede( -offset );
    return *this;
}
ConstVectorDataIterator &ConstVectorDataIterator::operator+=( int offset )
{
    if ( offset > 0 )
        advance( offset );
    if ( offset < 0 )
        recede( -offset );
    return *this;
}
VectorDataIterator &VectorDataIterator::operator-=( int offset )
{
    if ( offset > 0 )
        recede( offset );
    if ( offset < 0 )
        advance( -offset );
    return *this;
}
ConstVectorDataIterator &ConstVectorDataIterator::operator-=( int offset )
{
    if ( offset > 0 )
        recede( offset );
    if ( offset < 0 )
        advance( -offset );
    return *this;
}
VectorDataIterator VectorDataIterator::operator+( int offset )
{
    VectorDataIterator ans( *this );
    if ( offset > 0 )
        ans.advance( offset );
    if ( offset < 0 )
        ans.recede( -offset );
    return ans;
}
ConstVectorDataIterator ConstVectorDataIterator::operator+( int offset )
{
    ConstVectorDataIterator ans( *this );
    if ( offset > 0 )
        ans.advance( offset );
    if ( offset < 0 )
        ans.recede( -offset );
    return ans;
}
VectorDataIterator VectorDataIterator::operator-( int offset )
{
    VectorDataIterator ans( *this );
    if ( offset > 0 )
        ans.recede( offset );
    if ( offset < 0 )
        ans.advance( -offset );
    return ans;
}
ConstVectorDataIterator ConstVectorDataIterator::operator-( int offset )
{
    ConstVectorDataIterator ans( *this );
    if ( offset > 0 )
        ans.recede( offset );
    if ( offset < 0 )
        ans.advance( -offset );
    return ans;
}


/****************************************************************
* Difference Operators                                          *
****************************************************************/
int VectorDataIterator::operator-( const VectorDataIterator &rhs ) const
{
    return d_position - rhs.d_position;
}
int ConstVectorDataIterator::operator-( const ConstVectorDataIterator &rhs ) const
{
    return d_position - rhs.d_position;
}


/****************************************************************
* Equal Operators                                               *
****************************************************************/
bool VectorDataIterator::operator==( const VectorDataIterator &rhs ) const
{
    return d_Vec == rhs.d_Vec && d_position == rhs.d_position;
}
bool ConstVectorDataIterator::operator==( const ConstVectorDataIterator &rhs ) const
{
    return d_Vec == rhs.d_Vec && d_position == rhs.d_position;
}
bool VectorDataIterator::operator!=( const VectorDataIterator &rhs ) const
{
    return d_Vec != rhs.d_Vec || d_position != rhs.d_position;
}
bool ConstVectorDataIterator::operator!=( const ConstVectorDataIterator &rhs ) const
{
    return d_Vec != rhs.d_Vec || d_position != rhs.d_position;
}


/****************************************************************
* Assigment Operators                                           *
****************************************************************/
double &VectorDataIterator::operator[]( int i )
{
    VectorDataIterator tmp( *this ); // Create a temporary variable
    if ( i > 0 )
        tmp.advance( i );
    if ( i < 0 )
        tmp.recede( -i );
    return tmp.d_Block[tmp.d_CurOffset];
}
const double &ConstVectorDataIterator::operator[]( int i )
{
    ConstVectorDataIterator tmp( *this ); // Create a temporary variable
    if ( i > 0 )
        tmp.advance( i );
    if ( i < 0 )
        tmp.recede( -i );
    return tmp.d_Block[tmp.d_CurOffset];
}
double &VectorDataIterator::operator*() { return d_Block[d_CurOffset]; }
const double &ConstVectorDataIterator::operator*() { return d_Block[d_CurOffset]; }
}
}
