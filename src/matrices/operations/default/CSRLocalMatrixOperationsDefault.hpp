#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/CSRMatrixOperationsDefault.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Vector.h"

#include <algorithm>
#include <memory>

#include "ProfilerApp.h"

namespace AMP::LinearAlgebra {

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRLocalMatrixOperationsDefault<Policy, Allocator>::mult( const typename Policy::scalar_t *in,
                                                               std::shared_ptr<LocalMatrixData> A,
                                                               typename Policy::scalar_t *out )
{
    const auto nRows                   = static_cast<lidx_t>( csrData->numLocalRows() );
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
void CSRMatrixOperationsDefault<Policy, Allocator>::multTranspose(
    const typename Policy::scalar_t *in,
    std::shared_ptr<LocalMatrixData> A,
    std::vector<typename Policy::scalar_t> &vvals,
    std::vector<typename Policy::size_t> &vvals )
{

    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );

#error Below copies from diag only portion, update col count
    auto [nnz, cols, cols_loc, coeffs] = A->getDataFields();
    const auto num_unq                 = A->numLocalColumns();

    vvals.resize( num_unq, 0.0 );
    rcols.resize( num_unq );

    lidx_t offset = 0;
    for ( lidx_t row = 0; row < nRows; ++row ) {

        const auto ncols = nnz[row];
        const auto cloc  = &cols_loc[offset];
        const auto vloc  = &coeffs[offset];
        const auto val   = inDataBlock[row];

        for ( lidx_t j = 0; j < ncols; ++j ) {
            rcols[cloc[j]] = cols[offset + j];
            vvals[cloc[j]] += vloc[j] * val;
        }

        offset += ncols;
    }
}

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRMatrixOperationsDefault<Policy, Allocator>::scale( typename Policy::scalar_t alpha,
                                                           std::shared_ptr<LocalMatrixData> A )
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
void CSRMatrixOperationsDefault<Policy, Allocator>::matMultiply( std::shared_ptr<LocalMatrixData>,
                                                                 std::shared_ptr<LocalMatrixData>,
                                                                 std::shared_ptr<LocalMatrixData> )
{
    AMP_WARNING( "SpGEMM for CSRMatrixOperationsDefault not implemented" );
}

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRMatrixOperationsDefault<Policy, Allocator>::axpy( typename Policy::scalar_t alpha,
                                                          std::shared_ptr<LocalMatrixData> X,
                                                          std::shared_ptr<LocalMatrixData> Y )
{
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    const auto [nnz_x, cols_x, cols_loc_x, coeffs_x] = X->getDataFields();
    auto [nnz_y, cols_y, cols_loc_y, coeffs_y]       = Y->getDataFields();

    const auto tnnz = X->numberOfNonZeros();
    for ( lidx_t i = 0; i < tnnz; ++i ) {
        coeffs_y[i] += alpha * coeffs_x[i];
    }
}

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRMatrixOperationsDefault<Policy, Allocator>::setScalar( typename Policy::scalar_t alpha,
                                                               std::shared_ptr<LocalMatrixData> A )
{
    using scalar_t = typename Policy::scalar_t;

    auto [nnz, cols, cols_loc, coeffs] = A->getDataFields();

    const auto tnnz = A->numberOfNonZeros();

    std::fill( coeffs, coeffs + tnnz, alpha );
}

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRMatrixOperationsDefault<Policy, Allocator>::zero( std::shared_ptr<LocalMatrixData> A )
{
    using scalar_t = typename Policy::scalar_t;
    setScalar( static_cast<scalar_t>( 0.0 ), A );
}

template<typename Policy, class Allocator, class LocalMatrixData>
void CSRMatrixOperationsDefault<Policy, Allocator>::setDiagonal(
    const typename Policy::scalar_t *in, std::shared_ptr<LocalMatrixData> A )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

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
void CSRMatrixOperationsDefault<Policy, Allocator>::setIdentity(
    std::shared_ptr<LocalMatrixData> A )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto [nnz, cols, cols_loc, coeffs] = csrData->getDataFields();

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
void CSRMatrixOperationsDefault<Policy, Allocator>::extractDiagonal(
    std::shared_ptr<LocalMatrixData> A, typename Policy::scalar_t *buf )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto [nnz_d, cols_d, cols_loc_d, coeffs_d] = csrData->getDataFields();

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
void CSRMatrixOperationsDefault<Policy, Allocator>::LinfNorm(
    std::shared_ptr<LocalMatrixData> A, typename Policy::scalar_t *rowSums ) const
{
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto [nnz, cols, cols_loc, coeffs] = A->getCSRDiagData();
    auto rs                            = A->getDiagRowStarts();

    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );

    for ( lidx_t row = 0; row < nRows; ++row ) {
        auto nCols = nnz[row];
        auto start = rs[row];
        for ( lidx_t j = 0; j < nCols; ++j ) {
            rowSums[row] += std::abs( coeffs[start + j] );
        }
    }
}

} // namespace AMP::LinearAlgebra
