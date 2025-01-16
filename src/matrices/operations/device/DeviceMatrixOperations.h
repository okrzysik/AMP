#ifndef included_DeviceMatrixOperationsHelpers_H_
#define included_DeviceMatrixOperationsHelpers_H_

#include <cstddef>

namespace AMP {
namespace LinearAlgebra {


template<typename G, typename L, typename S>
struct DeviceMatrixOperations {
    static void mult( const L *row_starts,
                      const G *cols,
                      const S *coeffs,
                      const size_t N,
                      const S *in_h,
                      const size_t Ng,
                      S *out );
    static void mult( const L *row_starts,
                      const G *cols,
                      const S *coeffs,
                      const size_t N,
                      const S *in_h,
                      S *out );

    static void setScalar( const size_t N, S *coeffs, const S alpha );
    static void scale( const size_t N, S *coeffs, const S alpha );
    static void axpy( const size_t N, const S alpha, S *x, S *y );
    static void setDiagonal( const L *row_starts,
                             const G *cols,
                             S *coeffs,
                             const size_t N,
                             const size_t first_col,
                             const S *diag );
    static void extractDiagonal( const L *row_starts,
                                 const G *cols,
                                 const S *coeffs,
                                 const size_t N,
                                 const size_t first_col,
                                 S *diag );
    static void setIdentity(
        const L *row_starts, const G *cols, S *coeffs, const size_t N, const size_t first_col );
    static void LinfNorm( const size_t N, const S *x, const L *row_starts, S *row_sums );
};

} // namespace LinearAlgebra
} // namespace AMP
#endif
