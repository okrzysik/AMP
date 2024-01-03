#include "AMP/AMP_TPLs.h"

#if defined( AMP_USE_HYPRE )

    #include "AMP/matrices/CSRMatrix.hpp"
    #include "AMP/matrices/data/hypre/HypreCSRPolicy.h"

namespace AMP::LinearAlgebra {
template class CSRMatrix<HypreCSRPolicy>;
}

#endif
