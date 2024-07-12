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

#ifdef USE_HIP
    #include "AMP/utils/hip/HipAllocator.h"
namespace AMP::LinearAlgebra {
template class CSRMatrixOperationsDefault<CSRPolicy<size_t, int, double>, AMP::HipDevAllocator<int>>;
template class CSRMatrixData<CSRPolicy<size_t, int, double>, AMP::HipDevAllocator<int>>;
template class CSRMatrix<CSRPolicy<size_t, int, double>, AMP::HipDevAllocator<int>>;
template class CSRMatrixOperationsDefault<CSRPolicy<size_t, int, double>, AMP::HipManagedAllocator<int>>;
template class CSRMatrixData<CSRPolicy<size_t, int, double>, AMP::HipManagedAllocator<int>>;
template class CSRMatrix<CSRPolicy<size_t, int, double>, AMP::HipManagedAllocator<int>>;
} // namespace AMP::LinearAlgebra
#endif


#if defined( AMP_USE_HYPRE )

    #include "AMP/matrices/data/hypre/HypreCSRPolicy.h"

namespace AMP::LinearAlgebra {
template class CSRMatrixOperationsDefault<HypreCSRPolicy, std::allocator<int>>;
template class CSRMatrix<HypreCSRPolicy, std::allocator<int>>;
} // namespace AMP::LinearAlgebra

#endif
