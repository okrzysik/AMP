#ifndef included_AMP_MemoryPool_hpp
#define included_AMP_MemoryPool_hpp


#include "AMP/utils/threadpool/atomic_helpers.h"

#include <iostream>
#include <stdexcept>
#include <thread>


namespace AMP {


/******************************************************************
 * MemoryPool                                                      *
 ******************************************************************/
template<class TYPE, class INT_TYPE>
MemoryPool<TYPE, INT_TYPE>::MemoryPool( size_t size )
{
    static_assert( sizeof( TYPE ) >= sizeof( int ),
        "sizeof(TYPE) must be >= sizeof(int) to ensure proper operation" );
    static_assert( sizeof( TYPE ) >= sizeof( INT_TYPE ),
        "sizeof(TYPE) must be >= sizeof(INT_TYPE) to ensure proper operation" );
    d_objects = reinterpret_cast<TYPE *>( malloc( sizeof( TYPE ) * size ) );
    d_next    = 1;
    for ( size_t i = 0; i < size; i++ )
        reinterpret_cast<volatile INT_TYPE &>( d_objects[i] ) = i + 1;
    reinterpret_cast<volatile INT_TYPE &>( d_objects[size - 1] ) = -1;
}
template<class TYPE, class INT_TYPE>
MemoryPool<TYPE, INT_TYPE>::~MemoryPool()
{
    free( const_cast<TYPE *>( d_objects ) );
    d_objects = nullptr;
}
template<class TYPE, class INT_TYPE>
inline TYPE *MemoryPool<TYPE, INT_TYPE>::allocate()
{
    AtomicOperations::int32_atomic i = 0;
    while ( i == 0 )
        AtomicOperations::atomic_swap( &d_next, &i );
    TYPE *ptr = nullptr;
    if ( i != -1 ) {
        INT_TYPE j = reinterpret_cast<volatile INT_TYPE &>( d_objects[i - 1] );
        ptr        = const_cast<TYPE *>( &d_objects[i - 1] );
        new ( ptr ) TYPE();
        i = j + 1;
    }
    AtomicOperations::atomic_fetch_and_or( &d_next, i );
    return ptr;
}
template<class TYPE, class INT_TYPE>
inline void MemoryPool<TYPE, INT_TYPE>::free( TYPE *ptr )
{
    ptr->~TYPE();
    AtomicOperations::int32_atomic i = 0;
    while ( i == 0 )
        AtomicOperations::atomic_swap( &d_next, &i );
    reinterpret_cast<INT_TYPE &>( *ptr ) = i - 1;
    i                                    = ptr - d_objects + 1;
    AtomicOperations::atomic_fetch_and_or( &d_next, i );
}


} // namespace AMP

#endif
