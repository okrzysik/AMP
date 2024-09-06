#include "AMP/matrices/operations/CSRMatrixOperationsDefault.hpp"
#include "AMP/AMP_TPLs.h"
#include "AMP/matrices/CSRPolicy.h"
#include "AMP/utils/memory.h"

#define INSTANTIATE_FULL( policy, allocator ) \
    template class AMP::LinearAlgebra::CSRMatrixOperationsDefault<policy, allocator>;

// Check if device based allocators are needed
#ifdef USE_DEVICE
    #define INSTANTIATE_ALLOCS( policy )                       \
        INSTANTIATE_FULL( policy, AMP::HostAllocator<int> )    \
        INSTANTIATE_FULL( policy, AMP::ManagedAllocator<int> ) \
        INSTANTIATE_FULL( policy, AMP::DeviceAllocator<int> )
#else
    #define INSTANTIATE_ALLOCS( policy ) INSTANTIATE_FULL( policy, AMP::HostAllocator<int> )
#endif

// Check if hypre is present
using CSRPolicyDouble = AMP::LinearAlgebra::CSRPolicy<size_t, int, double>;
#if defined( AMP_USE_HYPRE )
    #include "AMP/matrices/data/hypre/HypreCSRPolicy.h"
    #define INSTANTIATE_POLICIES              \
        INSTANTIATE_ALLOCS( CSRPolicyDouble ) \
        INSTANTIATE_ALLOCS( AMP::LinearAlgebra::HypreCSRPolicy )
#else
    #define INSTANTIATE_POLICIES INSTANTIATE_ALLOCS( CSRPolicyDouble )
#endif

INSTANTIATE_POLICIES
