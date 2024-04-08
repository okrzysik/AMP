#include "AMP/matrices/CSRPolicy.h"
#include <cstddef>
#include <cstdint>

namespace AMP::LinearAlgebra {

template struct CSRPolicy<size_t, size_t, double>;
template struct CSRPolicy<size_t, size_t, float>;

template struct CSRPolicy<size_t, int32_t, double>;
template struct CSRPolicy<size_t, int32_t, float>;

} // namespace AMP::LinearAlgebra

#ifdef AMP_USE_HYPRE

extern "C" {
    #include "HYPRE.h"
    #include "HYPRE_IJ_mv.h"
    #include "HYPRE_utilities.h"
}

namespace AMP::LinearAlgebra {
template struct CSRPolicy<HYPRE_BigInt, HYPRE_Int, HYPRE_Real>;
}

#endif
