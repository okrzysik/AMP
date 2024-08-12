#include "AMP/matrices/CSRMatrix.hpp"
#include "AMP/AMP_TPLs.h"
#include "AMP/matrices/CSRPolicy.h"
#include "AMP/matrices/data/CSRMatrixData.hpp"
#include "AMP/matrices/operations/CSRMatrixOperationsDefault.hpp"
#include "AMP/utils/memory.h"

namespace AMP::LinearAlgebra {
template class CSRMatrixOperationsDefault<CSRPolicy<size_t, int, double>, std::allocator>;
template class CSRMatrixData<CSRPolicy<size_t, int, double>, std::allocator>;
template class CSRMatrix<CSRPolicy<size_t, int, double>, std::allocator>;
} // namespace AMP::LinearAlgebra

#if ( defined USE_DEVICE )
namespace AMP::LinearAlgebra {
template class CSRMatrixOperationsDefault<CSRPolicy<size_t, int, double>, AMP::DeviceAllocator>;
template class CSRMatrixData<CSRPolicy<size_t, int, double>, AMP::DeviceAllocator>;
template class CSRMatrix<CSRPolicy<size_t, int, double>, AMP::DeviceAllocator>;

template class CSRMatrixOperationsDefault<CSRPolicy<size_t, int, double>, AMP::ManagedAllocator>;
template class CSRMatrixData<CSRPolicy<size_t, int, double>, AMP::ManagedAllocator>;
template class CSRMatrix<CSRPolicy<size_t, int, double>, AMP::ManagedAllocator>;
} // namespace AMP::LinearAlgebra
#endif

#if defined( AMP_USE_HYPRE )

    #include "AMP/matrices/data/hypre/HypreCSRPolicy.h"

namespace AMP::LinearAlgebra {
template class CSRMatrixOperationsDefault<HypreCSRPolicy>;
template class CSRMatrix<HypreCSRPolicy>;
} // namespace AMP::LinearAlgebra

#endif
