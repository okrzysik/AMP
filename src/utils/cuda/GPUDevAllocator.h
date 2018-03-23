#ifndef included_AMP_GPUDevAllocator_H_
#define included_AMP_GPUDevAllocator_H_

template<typename T>
class GPUDevAllocator
{
public:
    T *allocate( size_t n );
    void deallocate( T *p, size_t n );
    void construct( T *p );
    void destroy( T *p );
};

#include "GPUDevAllocator.hpp"

#endif
