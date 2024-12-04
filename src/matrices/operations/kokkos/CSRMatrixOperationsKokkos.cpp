#include "AMP/matrices/operations/kokkos/CSRMatrixOperationsKokkos.hpp"
#include "AMP/AMP_TPLs.h"
#include "AMP/matrices/CSRPolicy.h"
#include "AMP/matrices/data/CSRLocalMatrixData.h"
#include "AMP/matrices/operations/kokkos/CSRLocalMatrixOperationsKokkos.hpp"
#include "AMP/utils/memory.h"

#if defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS )

    // This chain of macros is best read from bottom to top

    // Full instantiator macro is final one in the chain
    // this gives all template arguments
    #define INSTANTIATE_FULL( policy, allocator, execspace, viewspace )    \
        template class AMP::LinearAlgebra::CSRLocalMatrixOperationsKokkos< \
            policy,                                                        \
            allocator,                                                     \
            execspace,                                                     \
            viewspace,                                                     \
            AMP::LinearAlgebra::CSRLocalMatrixData<policy, allocator>>;    \
        template class AMP::LinearAlgebra::CSRMatrixOperationsKokkos<      \
            policy,                                                        \
            allocator,                                                     \
            execspace,                                                     \
            viewspace,                                                     \
            AMP::LinearAlgebra::CSRLocalMatrixData<policy, allocator>,     \
            AMP::LinearAlgebra::CSRLocalMatrixData<policy, allocator>>;

    // Execution spaces and view spaces are instatiated together and anchor
    // the chain. Valid space combinations are applied to the given
    // policy and allocator
    #ifdef USE_DEVICE
        #define INSTANTIATE_SPACES( policy, allocator )                                   \
            INSTANTIATE_FULL(                                                             \
                policy, allocator, Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace ) \
            INSTANTIATE_FULL(                                                             \
                policy, allocator, Kokkos::DefaultExecutionSpace, Kokkos::SharedSpace )   \
            INSTANTIATE_FULL( policy,                                                     \
                              allocator,                                                  \
                              Kokkos::DefaultExecutionSpace,                              \
                              Kokkos::DefaultExecutionSpace::memory_space )
    #else
        #define INSTANTIATE_SPACES( policy, allocator ) \
            INSTANTIATE_FULL(                           \
                policy, allocator, Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace )
    #endif

    // Allocator instatiator is responsible for adding all supported
    // allocator types for a given policy and forwarding along
    #ifdef USE_DEVICE
        #define INSTANTIATE_ALLOCS( policy )                          \
            INSTANTIATE_SPACES( policy, AMP::HostAllocator<void> )    \
            INSTANTIATE_SPACES( policy, AMP::ManagedAllocator<void> ) \
            INSTANTIATE_SPACES( policy, AMP::DeviceAllocator<void> )
    #else
        #define INSTANTIATE_ALLOCS( policy ) INSTANTIATE_SPACES( policy, AMP::HostAllocator<void> )
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

#endif
