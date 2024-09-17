#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/default/CSRLocalMatrixOperationsDefault.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Vector.h"

#include <algorithm>
#include <memory>

#include "ProfilerApp.h"

namespace AMP::LinearAlgebra {

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRLocalMatrixOperationsDefault<Policy, Allocator, LocalMatrixData>::mult(
    const typename Policy::scalar_t *in,
    std::shared_ptr<LocalMatrixData> A,
    typename Policy::scalar_t *out )
{
    using lidx_t = typename Policy::lidx_t;

    const auto nRows                   = static_cast<lidx_t>( A->numLocalRows() );
    auto [nnz, cols, cols_loc, coeffs] = A->getDataFields();
    lidx_t offset                      = 0;
    for ( lidx_t row = 0; row < nRows; ++row ) {
        const auto nCols = nnz[row];
        const auto cloc  = &cols_loc[offset];
        const auto vloc  = &coeffs[offset];
        for ( lidx_t c = 0; c < nCols; ++c ) {
            // Note: output is assumed to have useful values already
            out[row] += vloc[c] * in[cloc[c]];
        }

        offset += nCols;
    }
}

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRLocalMatrixOperationsDefault<Policy, Allocator, LocalMatrixData>::multTranspose(
    const typename Policy::scalar_t *in,
    std::shared_ptr<LocalMatrixData> A,
    std::vector<typename Policy::scalar_t> &vvals,
    std::vector<size_t> &rcols )
{
    using lidx_t = typename Policy::lidx_t;

    const auto nRows = static_cast<lidx_t>( A->numLocalRows() );

    auto [nnz, cols, cols_loc, coeffs] = A->getDataFields();
    A->getColumnMap( rcols );
    vvals.resize( rcols.size(), 0.0 );

    lidx_t offset = 0;
    for ( lidx_t row = 0; row < nRows; ++row ) {

        const auto ncols = nnz[row];
        const auto cloc  = &cols_loc[offset];
        const auto vloc  = &coeffs[offset];
        const auto val   = in[row];

        for ( lidx_t j = 0; j < ncols; ++j ) {
            vvals[cloc[j]] += vloc[j] * val;
        }

        offset += ncols;
    }
}

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRLocalMatrixOperationsDefault<Policy, Allocator, LocalMatrixData>::scale(
    typename Policy::scalar_t alpha, std::shared_ptr<LocalMatrixData> A )
{
    using scalar_t = typename Policy::scalar_t;

    auto [nnz, cols, cols_loc, coeffs] = A->getDataFields();

    const auto tnnz = A->numberOfNonZeros();

    std::transform( const_cast<scalar_t *>( coeffs ),
                    const_cast<scalar_t *>( coeffs ) + tnnz,
                    const_cast<scalar_t *>( coeffs ),
                    [=]( scalar_t val ) { return alpha * val; } );
}

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRLocalMatrixOperationsDefault<Policy, Allocator, LocalMatrixData>::matMultiply(
    std::shared_ptr<LocalMatrixData>,
    std::shared_ptr<LocalMatrixData>,
    std::shared_ptr<LocalMatrixData> )
{
    AMP_WARNING( "SpGEMM for CSRLocalMatrixOperationsDefault not implemented" );
}

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRLocalMatrixOperationsDefault<Policy, Allocator, LocalMatrixData>::axpy(
    typename Policy::scalar_t alpha,
    std::shared_ptr<LocalMatrixData> X,
    std::shared_ptr<LocalMatrixData> Y )
{
    using lidx_t = typename Policy::lidx_t;

    const auto [nnz_x, cols_x, cols_loc_x, coeffs_x] = X->getDataFields();
    auto [nnz_y, cols_y, cols_loc_y, coeffs_y]       = Y->getDataFields();

    const auto tnnz = X->numberOfNonZeros();
    for ( lidx_t i = 0; i < tnnz; ++i ) {
        coeffs_y[i] += alpha * coeffs_x[i];
    }
}

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRLocalMatrixOperationsDefault<Policy, Allocator, LocalMatrixData>::setScalar(
    typename Policy::scalar_t alpha, std::shared_ptr<LocalMatrixData> A )
{
    auto [nnz, cols, cols_loc, coeffs] = A->getDataFields();
    const auto tnnz                    = A->numberOfNonZeros();
    std::fill( coeffs, coeffs + tnnz, alpha );
}

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRLocalMatrixOperationsDefault<Policy, Allocator, LocalMatrixData>::zero(
    std::shared_ptr<LocalMatrixData> A )
{
    using scalar_t = typename Policy::scalar_t;
    setScalar( static_cast<scalar_t>( 0.0 ), A );
}

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRLocalMatrixOperationsDefault<Policy, Allocator, LocalMatrixData>::setDiagonal(
    const typename Policy::scalar_t *in, std::shared_ptr<LocalMatrixData> A )
{
    using lidx_t = typename Policy::lidx_t;
    using gidx_t = typename Policy::gidx_t;

    const auto nRows    = static_cast<lidx_t>( A->numLocalRows() );
    const auto beginRow = A->beginRow();

    auto [nnz, cols, cols_loc, coeffs] = A->getDataFields();

    lidx_t offset = 0;
    for ( lidx_t row = 0; row < nRows; ++row ) {
        const auto ncols = nnz[row];
        for ( lidx_t icol = 0; icol < ncols; ++icol ) {
            if ( cols[offset + icol] == static_cast<gidx_t>( beginRow + row ) ) {
                coeffs[offset + icol] = in[row];
                break;
            }
        }
        offset += nnz[row];
    }
}

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRLocalMatrixOperationsDefault<Policy, Allocator, LocalMatrixData>::setIdentity(
    std::shared_ptr<LocalMatrixData> A )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto [nnz, cols, cols_loc, coeffs] = A->getDataFields();

    const auto nRows    = static_cast<lidx_t>( A->numLocalRows() );
    const auto beginRow = A->beginRow();

    lidx_t offset = 0;
    for ( lidx_t row = 0; row < nRows; ++row ) {
        const auto ncols = nnz[row];
        for ( lidx_t icol = 0; icol < ncols; ++icol ) {
            if ( cols[offset + icol] == static_cast<gidx_t>( beginRow + row ) ) {
                coeffs[offset + icol] = static_cast<scalar_t>( 1.0 );
                break;
            }
        }
        offset += nnz[row];
    }
}

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRLocalMatrixOperationsDefault<Policy, Allocator, LocalMatrixData>::extractDiagonal(
    std::shared_ptr<LocalMatrixData> A, typename Policy::scalar_t *buf )
{
    using lidx_t = typename Policy::lidx_t;
    using gidx_t = typename Policy::gidx_t;

    auto [nnz, cols, cols_loc, coeffs] = A->getDataFields();

    const auto nRows    = static_cast<lidx_t>( A->numLocalRows() );
    const auto beginRow = A->beginRow();

    lidx_t offset = 0;
    for ( lidx_t row = 0; row < nRows; ++row ) {
        const auto ncols = nnz[row];
        for ( lidx_t icol = 0; icol < ncols; ++icol ) {
            if ( cols[offset + icol] == static_cast<gidx_t>( beginRow + row ) ) {
                buf[row] = coeffs[offset + icol];
                break;
            }
        }
        offset += nnz[row];
    }
}

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRLocalMatrixOperationsDefault<Policy, Allocator, LocalMatrixData>::LinfNorm(
    std::shared_ptr<LocalMatrixData> A, typename Policy::scalar_t *rowSums ) const
{
    using lidx_t = typename Policy::lidx_t;

    auto [nnz, cols, cols_loc, coeffs] = A->getDataFields();
    auto rs                            = A->getRowStarts();

    const auto nRows = static_cast<lidx_t>( A->numLocalRows() );

    for ( lidx_t row = 0; row < nRows; ++row ) {
        auto nCols = nnz[row];
        auto start = rs[row];
        for ( lidx_t j = 0; j < nCols; ++j ) {
            rowSums[row] += std::abs( coeffs[start + j] );
        }
    }
}

} // namespace AMP::LinearAlgebra
