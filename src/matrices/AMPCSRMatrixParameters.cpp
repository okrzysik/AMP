#include "AMP/matrices/AMPCSRMatrixParameters.hpp"
#include "AMP/AMP_TPLs.h"
#include "AMP/matrices/CSRPolicy.h"

#define INSTANTIATE_FULL( policy ) \
    template class AMP::LinearAlgebra::AMPCSRMatrixParameters<policy>;

// Check if hypre is present
using CSRPolicyDouble = AMP::LinearAlgebra::CSRPolicy<size_t, int, double>;
#if defined( AMP_USE_HYPRE )
    #include "AMP/matrices/data/hypre/HypreCSRPolicy.h"
    #define INSTANTIATE_POLICIES            \
        INSTANTIATE_FULL( CSRPolicyDouble ) \
        INSTANTIATE_FULL( AMP::LinearAlgebra::HypreCSRPolicy )
#else
    #define INSTANTIATE_POLICIES INSTANTIATE_FULL( CSRPolicyDouble )
#endif

INSTANTIATE_POLICIES
