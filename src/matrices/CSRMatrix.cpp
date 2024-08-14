#include "AMP/matrices/CSRMatrix.hpp"
#include "AMP/AMP_TPLs.h"
#include "AMP/matrices/CSRPolicy.h"
#include "AMP/matrices/data/CSRMatrixData.hpp"
#include "AMP/matrices/operations/CSRMatrixOperationsDefault.hpp"
#include "AMP/matrices/operations/CSRMatrixOperationsKokkos.hpp"
#include "AMP/utils/memory.h"

namespace AMP::LinearAlgebra {
template class CSRMatrixOperationsDefault<CSRPolicy<size_t, int, double>, AMP::HostAllocator<int>>;
template class CSRMatrixData<CSRPolicy<size_t, int, double>, AMP::HostAllocator<int>>;
template class CSRMatrix<CSRPolicy<size_t, int, double>, AMP::HostAllocator<int>>;
} // namespace AMP::LinearAlgebra

#if ( defined USE_DEVICE )
namespace AMP::LinearAlgebra {
template class CSRMatrixOperationsDefault<CSRPolicy<size_t, int, double>,
                                          AMP::DeviceAllocator<int>>;
template class CSRMatrixData<CSRPolicy<size_t, int, double>, AMP::DeviceAllocator<int>>;
template class CSRMatrix<CSRPolicy<size_t, int, double>, AMP::DeviceAllocator<int>>;

template class CSRMatrixOperationsDefault<CSRPolicy<size_t, int, double>,
                                          AMP::ManagedAllocator<int>>;
template class CSRMatrixData<CSRPolicy<size_t, int, double>, AMP::ManagedAllocator<int>>;
template class CSRMatrix<CSRPolicy<size_t, int, double>, AMP::ManagedAllocator<int>>;
} // namespace AMP::LinearAlgebra
#endif

#if defined( AMP_USE_HYPRE )

    #include "AMP/matrices/data/hypre/HypreCSRPolicy.h"

namespace AMP::LinearAlgebra {
template class CSRMatrixOperationsDefault<HypreCSRPolicy, AMP::HostAllocator<int>>;
template class CSRMatrix<HypreCSRPolicy, AMP::HostAllocator<int>>;


    #if ( defined USE_DEVICE )
template class CSRMatrixOperationsDefault<HypreCSRPolicy, AMP::DeviceAllocator<int>>;
template class CSRMatrixData<HypreCSRPolicy, AMP::DeviceAllocator<int>>;
template class CSRMatrix<HypreCSRPolicy, AMP::DeviceAllocator<int>>;

template class CSRMatrixOperationsDefault<HypreCSRPolicy, AMP::ManagedAllocator<int>>;
template class CSRMatrixData<HypreCSRPolicy, AMP::ManagedAllocator<int>>;
template class CSRMatrix<HypreCSRPolicy, AMP::ManagedAllocator<int>>;
    #endif
} // namespace AMP::LinearAlgebra

#endif
