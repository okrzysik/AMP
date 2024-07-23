#include "AMP/matrices/CSRMatrix.hpp"
#include "AMP/AMP_TPLs.h"
#include "AMP/matrices/CSRPolicy.h"
#include "AMP/matrices/data/CSRMatrixData.hpp"
#include "AMP/matrices/operations/CSRMatrixOperationsDefault.hpp"

namespace AMP::LinearAlgebra {
template class CSRMatrixOperationsDefault<CSRPolicy<size_t, int, double>, std::allocator<int>>;
template class CSRMatrixData<CSRPolicy<size_t, int, double>, std::allocator<int>>;
template class CSRMatrix<CSRPolicy<size_t, int, double>, std::allocator<int>>;
} // namespace AMP::LinearAlgebra

#include "AMP/utils/memory.h"
#ifdef USE_HIP
namespace AMP::LinearAlgebra {
template class CSRMatrixOperationsDefault<CSRPolicy<size_t, int, double>, AMP::DeviceAllocator<int>>;
template class CSRMatrixData<CSRPolicy<size_t, int, double>, AMP::DeviceAllocator<int>>;
template class CSRMatrix<CSRPolicy<size_t, int, double>, AMP::DeviceAllocator<int>>;
template class CSRMatrixOperationsDefault<CSRPolicy<size_t, int, double>, AMP::ManagedAllocator<int>>;
template class CSRMatrixData<CSRPolicy<size_t, int, double>, AMP::ManagedAllocator<int>>;
template class CSRMatrix<CSRPolicy<size_t, int, double>, AMP::ManagedAllocator<int>>;
} // namespace AMP::LinearAlgebra
#endif
#ifdef USE_CUDA
namespace AMP::LinearAlgebra {
template class CSRMatrixOperationsDefault<CSRPolicy<size_t, int, double>, AMP::DeviceAllocator<int>>;
template class CSRMatrixData<CSRPolicy<size_t, int, double>, AMP::DeviceAllocator<int>>;
template class CSRMatrix<CSRPolicy<size_t, int, double>, AMP::DeviceAllocator<int>>;
template class CSRMatrixOperationsDefault<CSRPolicy<size_t, int, double>, AMP::ManagedAllocator<int>>;
template class CSRMatrixData<CSRPolicy<size_t, int, double>, AMP::ManagedAllocator<int>>;
template class CSRMatrix<CSRPolicy<size_t, int, double>, AMP::ManagedAllocator<int>>;
} // namespace AMP::LinearAlgebra
#endif


#if defined( AMP_USE_HYPRE )

    #include "AMP/matrices/data/hypre/HypreCSRPolicy.h"

namespace AMP::LinearAlgebra {
template class CSRMatrixOperationsDefault<HypreCSRPolicy, std::allocator<int>>;
template class CSRMatrix<HypreCSRPolicy, std::allocator<int>>;
} // namespace AMP::LinearAlgebra

#endif
