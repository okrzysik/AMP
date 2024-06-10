#include "AMP/matrices/CSRMatrix.hpp"
#include "AMP/AMP_TPLs.h"
#include "AMP/matrices/CSRPolicy.h"
#include "AMP/matrices/data/CSRMatrixData.hpp"
#include "AMP/matrices/operations/CSRMatrixOperationsDefault.hpp"
#include "AMP/matrices/operations/CSRMatrixOperationsKokkos.hpp"

namespace AMP::LinearAlgebra {
template class CSRMatrixOperationsDefault<CSRPolicy<size_t, int, double>>;
template class CSRMatrixData<CSRPolicy<size_t, int, double>>;
template class CSRMatrix<CSRPolicy<size_t, int, double>>;
} // namespace AMP::LinearAlgebra

#if defined( AMP_USE_HYPRE )

    #include "AMP/matrices/data/hypre/HypreCSRPolicy.h"

namespace AMP::LinearAlgebra {
template class CSRMatrixOperationsDefault<HypreCSRPolicy>;
template class CSRMatrix<HypreCSRPolicy>;
} // namespace AMP::LinearAlgebra

#endif
