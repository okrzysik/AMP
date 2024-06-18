#include "AMP/utils/StackAllocator.h"

#include <cstring>
#include <iostream>
#include <utility>


/********************************************************************
 * Constructor/destructor                                            *
 ********************************************************************/
StackAllocator::StackAllocator( size_t bytes,
                                size_t blockSize,
                                std::function<void *( size_t )> allocator,
                                std::function<void( void * )> deallocator )
    : d_blockSize( 0 ),
      d_available( 0 ),
      d_N( 0 ),
      d_capacity( 1024 ),
      d_memory( nullptr ),
      d_ptr( nullptr ),
      d_allocator( std::move( allocator ) ),
      d_deallocator( std::move( deallocator ) )
{
    // Check that we have < 2^32 blocks
    if ( bytes / blockSize > 0xFFFFFFF0 )
        throw std::logic_error( "Maximum number of blocks exceeded, increase blockSize" );
    // Get the block size as 2^n
    d_blockSize = 0;
    while ( ( (size_t) 0x01 << d_blockSize ) < blockSize )
        d_blockSize++;
    if ( ( (size_t) 0x01 << d_blockSize ) != blockSize )
        throw std::logic_error( "blockSize must be a power of 2" );
    // Allocate the data
    d_memory = d_allocator( bytes );
    d_ptr    = new void *[d_capacity];
    for ( size_t i = 0; i < d_capacity; i++ )
        d_ptr[i] = nullptr;
    // Align the memory pointer to the block size
    d_ptr[0] = (void *) ( blockSize * ( ( (size_t) d_memory + blockSize - 1 ) / blockSize ) );
    // Get the number of blocks availible to allocate
    d_available = ( bytes - ( (size_t) d_ptr[0] - (size_t) d_memory ) ) / blockSize;
}
StackAllocator::StackAllocator()
    : d_blockSize( 0 ),
      d_available( 0 ),
      d_N( 0 ),
      d_capacity( 1024 ),
      d_memory( nullptr ),
      d_ptr( nullptr )
{
}
StackAllocator::StackAllocator( StackAllocator &&rhs )
    : d_blockSize( rhs.d_blockSize ),
      d_available( rhs.d_available ),
      d_N( rhs.d_N ),
      d_capacity( rhs.d_capacity ),
      d_memory( rhs.d_memory ),
      d_ptr( rhs.d_ptr ),
      d_allocator( std::move( rhs.d_allocator ) ),
      d_deallocator( std::move( rhs.d_deallocator ) )
{
    rhs.d_memory = nullptr;
}
StackAllocator &StackAllocator::operator=( StackAllocator &&rhs )
{
    if ( this == &rhs )
        return *this;
    std::swap( d_blockSize, rhs.d_blockSize );
    std::swap( d_available, rhs.d_available );
    std::swap( d_N, rhs.d_N );
    std::swap( d_capacity, rhs.d_capacity );
    std::swap( d_memory, rhs.d_memory );
    std::swap( d_ptr, rhs.d_ptr );
    std::swap( d_allocator, rhs.d_allocator );
    std::swap( d_deallocator, rhs.d_deallocator );
    return *this;
}
StackAllocator::~StackAllocator()
{
    if ( d_memory != nullptr ) {
        if ( d_N != 0 )
            std::cerr << "Some memory was not free'd before destroying StackAllocator\n";
        d_deallocator( d_memory );
        delete[] d_ptr;
    }
    d_memory = nullptr;
}


/********************************************************************
 * Allocate memory                                                   *
 ********************************************************************/
void *StackAllocator::allocate( size_t bytes )
{
    // Get the number of blocks needed
    size_t s = ( (size_t) 0x01 << d_blockSize );
    size_t n = ( bytes + s - 1 ) >> d_blockSize;
    if ( n == 0 || n > d_available )
        return nullptr;
    // Check that we can perform an allocation
    if ( d_N + 2 > d_capacity ) {
        void **tmp = d_ptr;
        d_capacity *= 2;
        d_ptr = new void *[d_capacity];
        memcpy( d_ptr, tmp, d_N + 1 );
        delete[] tmp;
    }
    // Perform the allocation
    d_ptr[d_N + 1] = (void *) ( (size_t) d_ptr[d_N] + ( n << d_blockSize ) );
    d_available -= n;
    d_N++;
    return d_ptr[d_N];
}


/********************************************************************
 * Free memory                                                       *
 ********************************************************************/
void StackAllocator::deallocate( void *ptr, std::size_t )
{
    if ( ptr == nullptr )
        return;
    if ( ptr != d_ptr[d_N] )
        throw std::bad_alloc();
    // Remove the allocation
    size_t n = ( (size_t) d_ptr[d_N] - (size_t) d_ptr[d_N - 1] ) >> d_blockSize;
    d_available += n;
    d_N--;
}
