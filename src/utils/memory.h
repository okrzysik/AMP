#ifndef included_AMP_MEMORY
#define included_AMP_MEMORY

#ifdef USE_CUDA
    #include "AMP/utils/cuda/CudaAllocator.h"
#endif
#ifdef USE_HIP
    #include "AMP/utils/hip/HipAllocator.h"
#endif

namespace AMP{
#ifdef USE_CUDA
    template<typename TYPE>
    using DeviceAllocator = AMP::CudaDevAllocator<TYPE>;
    template<typename TYPE>
    using ManagedAllocator = AMP::CudaManagedAllocator<TYPE>;
#endif
#ifdef USE_HIP
    template<typename TYPE>
    using DeviceAllocator = AMP::HipDevAllocator<TYPE>;
    template<typename TYPE>
    using ManagedAllocator = AMP::HipManagedAllocator<TYPE>;
#endif

template<typename TYPE>
using HostAllocator = std::allocator<TYPE>;
}
#endif
