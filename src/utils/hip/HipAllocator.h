#ifndef included_AMP_GPUDevAllocator
#define included_AMP_GPUDevAllocator

#include <hip/hip_runtime.h>
#include "AMP/utils/hip/helper_hip.h"

namespace AMP {

/**
 * \class  HipDevAllocator
 * @brief  Allocator based on hipMalloc
 */
template<typename T>
class HipDevAllocator
{
public:
    using value_type = T;

    T *allocate( size_t n )
    {
        T *ptr;
        auto err = hipMalloc( &ptr, n * sizeof( T ) );
        checkHipErrors( err );
        return ptr;
    }

    void deallocate( T *p, size_t )
    {
        auto err = hipFree( p );
        checkHipErrors( err );
    }
};


/**
 * \class  HipManagedAllocator
 * @brief  Allocator based on hipMallocManaged
 */
template<typename T>
class HipManagedAllocator
{
public:
    using value_type = T;

    T *allocate( size_t n )
    {
        T *ptr;
        auto err = hipMallocManaged( &ptr, n * sizeof( T ), hipMemAttachGlobal );
        checkHipErrors( err );
        return ptr;
    }

    void deallocate( T *p, size_t )
    {
        auto err = hipFree( p );
        checkHipErrors( err );
    }
};


/**
 * \class  HipHostAllocator
 * @brief  Allocator based on hipMallocHost
 */
template<typename T>
class HipHostAllocator
{
public:
    using value_type = T;

    T *allocate( size_t n )
    {
        T *ptr;
        auto err = hipHostMalloc( &ptr, n * sizeof( T ) );
        checkHipErrors( err );
        return ptr;
    }

    void deallocate( T *p, size_t )
    {
        auto err = hipHostFree( p );
        checkHipErrors( err );
    }
};

} // namespace AMP


#endif
