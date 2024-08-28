#ifndef included_AMP_MEMORY
#define included_AMP_MEMORY

#ifdef USE_CUDA
    #include "AMP/utils/cuda/CudaAllocator.h"
#endif
#ifdef USE_HIP
    #include "AMP/utils/hip/HipAllocator.h"
#endif

namespace AMP {
// managed allocators
#ifdef USE_CUDA
template<typename TYPE>
using ManagedAllocator = AMP::CudaManagedAllocator<TYPE>;
#elif defined( USE_HIP )
template<typename TYPE>
using ManagedAllocator = AMP::HipManagedAllocator<TYPE>;
#endif

// device allocators
#ifdef USE_CUDA
template<typename TYPE>
using DeviceAllocator = AMP::CudaDevAllocator<TYPE>;
#elif defined( USE_HIP )
template<typename TYPE>
using DeviceAllocator = AMP::HipDevAllocator<TYPE>;
#endif

// host allocators
#ifdef USE_CUDA
template<typename TYPE>
using HostAllocator = AMP::CudaHostAllocator<TYPE>;
#elif defined( USE_HIP )
template<typename TYPE>
using HostAllocator = AMP::HipHostAllocator<TYPE>;
#else // no device
template<typename TYPE>
using HostAllocator = std::allocator<TYPE>;
#endif

} // namespace AMP
#endif
