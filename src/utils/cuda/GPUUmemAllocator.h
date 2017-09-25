#ifndef included_GPUUmemAllocator_H_
#define included_GPUUmemAllocator_H_


template<typename T>
class GPUUmemAllocator
{
public:
    T *allocate( size_t n );
    void deallocate( T *p, size_t n );
    void construct( T *p );
    void destroy( T *p );
};

#include "GPUUmemAllocator.hpp"

#endif
