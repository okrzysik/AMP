#ifndef included_AMP_DevCopyCast_HPP_
#define included_AMP_DevCopyCast_HPP_

#include "AMP/utils/memory.h"
#include <memory>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/host_vector.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/mr/allocator.h>

#include <iostream>

namespace AMP::Utilities {

template<typename T1, typename T2>
struct copyCast_<T1, T2, AMP::Utilities::MemoryType::host> {
    void static apply( size_t len, const T1 *vec_in, T2 *vec_out )
    {
        auto lambda = [] __host__ __device__( T1 x ) { return static_cast<T2>( x ); };
        thrust::transform( thrust::host, vec_in, vec_in + len, vec_out, lambda );
    }
};

template<typename T1, typename T2>
struct copyCast_<T1, T2, AMP::Utilities::MemoryType::unregistered> {
    void static apply( size_t len, const T1 *vec_in, T2 *vec_out )
    {
        auto lambda = [] __host__ __device__( T1 x ) { return static_cast<T2>( x ); };
        thrust::transform( thrust::host, vec_in, vec_in + len, vec_out, lambda );
    }
};

template<typename T1, typename T2>
struct copyCast_<T1, T2, AMP::Utilities::MemoryType::managed> {
    void static apply( size_t len, const T1 *vec_in, T2 *vec_out )
    {
        auto lambda = [] __host__ __device__( T1 x ) { return static_cast<T2>( x ); };
        thrust::transform( thrust::device, vec_in, vec_in + len, vec_out, lambda );
    }
};

template<typename T1, typename T2>
struct copyCast_<T1, T2, AMP::Utilities::MemoryType::device> {
    void static apply( size_t len, const T1 *vec_in, T2 *vec_out )
    {
        auto lambda = [] __host__ __device__( T1 x ) { return static_cast<T2>( x ); };
        thrust::transform( thrust::device, vec_in, vec_in + len, vec_out, lambda );
    }
};

} // namespace AMP::Utilities

// #include "CopyCast.hpp"

#endif
