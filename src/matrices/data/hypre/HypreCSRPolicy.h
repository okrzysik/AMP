#ifndef included_HypreCSRPolicy_h
#define included_HypreCSRPolicy_h

extern "C" {
#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_utilities.h"
}

namespace AMP::LinearAlgebra {

struct HypreCSRPolicy {
    using gidx_t   = HYPRE_BigInt;
    using lidx_t   = HYPRE_Int;
    using scalar_t = HYPRE_Real;
};

} // namespace AMP::LinearAlgebra

#endif
