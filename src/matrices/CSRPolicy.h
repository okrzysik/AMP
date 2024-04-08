#ifndef included_CSRPolicy_H_
#define included_CSRPolicy_H_

namespace AMP::LinearAlgebra {

// G refers to the type for global indices
// L refers to the type for local indices
// S refers to the type for scalar values
template<typename G, typename L, typename S>
struct CSRPolicy {
    using gidx_t   = G;
    using lidx_t   = L;
    using scalar_t = S;
};
} // namespace AMP::LinearAlgebra

#endif
