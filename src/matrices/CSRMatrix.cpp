#include "AMP/AMP_TPLs.h"

#if defined( AMP_USE_HYPRE )

    #include "AMP/matrices/CSRMatrix.hpp"
    #include "AMP/matrices/data/hypre/HypreCSRPolicy.h"
    #include "AMP/matrices/operations/CSRMatrixOperationsDefault.hpp"

namespace AMP::LinearAlgebra {
template class CSRMatrixOperationsDefault<HypreCSRPolicy>;
template class CSRMatrix<HypreCSRPolicy>;
} // namespace AMP::LinearAlgebra

#endif
