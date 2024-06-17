#include "AMP/utils/BuddyAllocator.h"
#include "AMP/utils/StackAllocator.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/typeid.h"

#include <array>
#include <cmath>


// Dummy class for malloc
template<class TYPE>
class MallocAllocator
{
public:
    MallocAllocator(){};
    ~MallocAllocator(){};
    inline TYPE *allocate( size_t N ) { return (TYPE *) ::malloc( N * sizeof( TYPE ) ); }
    inline void deallocate( TYPE *ptr, size_t ) { ::free( ptr ); }
};


// Test if we are LIFO
template<typename T>
struct has_isLIFO {
private:
    typedef std::true_type yes;
    typedef std::false_type no;
    template<typename U>
    static auto test( int ) -> decltype( std::declval<U>().isLIFO() == 1, yes() );
    template<typename>
    static no test( ... );

public:
    static constexpr bool value = std::is_same<decltype( test<T>( 0 ) ), yes>::value;
};
template<class Allocator>
constexpr bool isLIFO()
{
    if constexpr ( has_isLIFO<Allocator>::value ) {
        return Allocator::isLIFO();
    } else {
        return false;
    }
}


// Test a memory allocator
template<class MEMORY, class TYPE>
void testAllocator( MEMORY &alloc, TYPE &size, AMP::UnitTest &ut )
{
    std::string name = AMP::getTypeID<MEMORY>().name;

    // Run some simple checks
    auto x    = alloc.allocate( 1 );
    auto y    = alloc.allocate( 1 );
    auto z    = alloc.allocate( 1 );
    bool pass = x != nullptr && y != nullptr && z != nullptr;
    alloc.deallocate( z, 1 );
    alloc.deallocate( y, 1 );
    alloc.deallocate( x, 1 );
    if ( !isLIFO<MEMORY>() ) {
        x = alloc.allocate( 1 );
        y = alloc.allocate( 1 );
        z = alloc.allocate( 1 );
        alloc.deallocate( x, 1 );
        alloc.deallocate( y, 1 );
        alloc.deallocate( z, 1 );
    }

    // Test performance
    printf( "%28s", name.data() );
    for ( size_t i = 0; i < size.size(); i++ ) {
        int N_it = 30000;
        x        = alloc.allocate( size[i] );
        alloc.deallocate( x, size[i] );
        auto t1 = std::chrono::high_resolution_clock::now();
        for ( int j = 0; j < N_it; j++ ) {
            y = alloc.allocate( size[i] );
            AMP::Utilities::nullUse( y );
            pass = pass && y != nullptr;
            alloc.deallocate( y, size[i] );
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        int ns  = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count() / N_it;
        printf( "%8i", ns );
    }
    printf( "\n" );

    // Finished
    if ( pass )
        ut.passes( name );
    else
        ut.failure( name );
}


// Main
int main( int, char *[] )
{
    AMP::UnitTest ut;

    // Total number of bytes and block size to allocate
    size_t block = 8;
    size_t bytes = 0x8000000; // 128 MB

    // Create output table
    std::array<size_t, 5> size = { 8, 1024, 0x100000, 0x2000000, bytes };
    printf( "%28s  8 bytes   1 kB    1 MB   32 MB   %i MB\n", "", (int) bytes / 0x100000 );


    // Test BuddyAllocator
    BuddyAllocator buddy( bytes, block );
    testAllocator( buddy, size, ut );
    buddy = BuddyAllocator();

    // Test StackAllocator
    StackAllocator stack( bytes, block );
    testAllocator( stack, size, ut );
    stack = StackAllocator();

    // Test malloc wrapper
    MallocAllocator<std::byte> malloc;
    testAllocator( malloc, size, ut );

    // Test std::allocator
    std::allocator<std::byte> standard;
    testAllocator( standard, size, ut );

    // Finished
    int N_errors = ut.NumFailGlobal();
    ut.report();
    ut.reset();
    if ( N_errors == 0 )
        std::cout << "All tests passed\n";
    return N_errors;
}
