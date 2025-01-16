#include "AMP/matrices/operations/device/CSRMatrixOperationsDevice.hpp"
#include "AMP/AMP_TPLs.h"
#include "AMP/matrices/CSRPolicy.h"
#include "AMP/matrices/data/CSRLocalMatrixData.h"
#include "AMP/matrices/operations/device/CSRLocalMatrixOperationsDevice.hpp"
#include "AMP/utils/memory.h"

// This chain of macros is best read from bottom to top

// Full instantiator macro is final one in the chain
// this gives all template arguments
#define INSTANTIATE_FULL( policy, allocator )                          \
    template class AMP::LinearAlgebra::CSRLocalMatrixOperationsDevice< \
        policy,                                                        \
        allocator,                                                     \
        AMP::LinearAlgebra::CSRLocalMatrixData<policy, allocator>>;    \
    template class AMP::LinearAlgebra::CSRMatrixOperationsDevice<      \
        policy,                                                        \
        allocator,                                                     \
        AMP::LinearAlgebra::CSRLocalMatrixData<policy, allocator>,     \
        AMP::LinearAlgebra::CSRLocalMatrixData<policy, allocator>>;

// Allocator instatiator is responsible for adding all supported
// allocator types for a given policy and forwarding along
#ifdef USE_DEVICE
    #define INSTANTIATE_ALLOCS( policy )                        \
        INSTANTIATE_FULL( policy, AMP::HostAllocator<void> )    \
        INSTANTIATE_FULL( policy, AMP::ManagedAllocator<void> ) \
        INSTANTIATE_FULL( policy, AMP::DeviceAllocator<void> )
#else
    #define INSTANTIATE_ALLOCS( policy ) INSTANTIATE_FULL( policy, AMP::HostAllocator<void> )
#endif

// Policy instantiator starts the chain and forwards
// all known policies to the next instantiator
// using CSRPolicyFloat  = AMP::LinearAlgebra::CSRPolicy<size_t, int, float>;
using CSRPolicyDouble = AMP::LinearAlgebra::CSRPolicy<size_t, int, double>;
#ifdef AMP_USE_HYPRE
    #include "AMP/matrices/data/hypre/HypreCSRPolicy.h"
    #define INSTANTIATE_POLICIES              \
        INSTANTIATE_ALLOCS( CSRPolicyDouble ) \
        INSTANTIATE_ALLOCS( AMP::LinearAlgebra::HypreCSRPolicy )
#else
    #define INSTANTIATE_POLICIES INSTANTIATE_ALLOCS( CSRPolicyDouble )
#endif

INSTANTIATE_POLICIES
