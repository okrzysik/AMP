#ifndef included_HypreCSRPolicy_h
#define included_HypreCSRPolicy_h

#include "AMP/matrices/CSRPolicy.h"

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_utilities.h"

namespace AMP::LinearAlgebra {
using HypreCSRPolicy = CSRPolicy<HYPRE_BigInt, HYPRE_Int, HYPRE_Real>;

} // namespace AMP::LinearAlgebra

#endif
