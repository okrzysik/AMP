#ifndef included_StackAllocator
#define included_StackAllocator


#include <cstdint>
#include <functional>
#include <stdlib.h>


/** \class StackAllocator
 *
 * This class provides basic routines to allocate/deallocate memory.
 * This allocator works like a stack and requires that items are free'd
 * in the reverse order in which they are created.
 * This class is not thread-safe.
 */
class StackAllocator
{
public:
    //! Default constructor
    explicit StackAllocator( size_t bytes,
                             size_t blockSize                          = 1024,
                             std::function<void *( size_t )> allocator = ::malloc,
                             std::function<void( void * )> deallocator = ::free );

    //! Empty constructor
    StackAllocator();

    //! Destructor
    ~StackAllocator();

    // Copy/assignment constructors
    StackAllocator( const StackAllocator & )            = delete;
    StackAllocator &operator=( const StackAllocator & ) = delete;
    StackAllocator( StackAllocator && );
    StackAllocator &operator=( StackAllocator && );

    //! Allocate memory
    void *allocate( size_t bytes );

    //! Deallocate memory
    void deallocate( void *p, size_t bytes );

    //! Check if the allocator is Last-In-First-Out (LIFO)
    static constexpr bool isLIFO() { return true; }

private:
    // Member data
    uint8_t d_blockSize;  // Blocksize stored as 2^n
    uint32_t d_available; // Number of blocks availible
    uint32_t d_N;         // Number of current allocations
    uint32_t d_capacity;  // Number of available allocations
    void *d_memory;       // Raw memory pointer
    void **d_ptr;         // Array of pointers allocated
    std::function<void *( size_t )> d_allocator;
    std::function<void( void * )> d_deallocator;
};


#endif
