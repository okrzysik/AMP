#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/CSRMatrixOperationsDefault.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>

namespace AMP::LinearAlgebra {

template<typename Policy>
static CSRMatrixData<Policy> const *getCSRMatrixData( MatrixData const &A )
{
    auto ptr = dynamic_cast<CSRMatrixData<Policy> const *>( &A );
    AMP_INSIST( ptr, "dynamic cast from const MatrixData to const CSRMatrixData failed" );
    return ptr;
}

template<typename Policy>
static CSRMatrixData<Policy> *getCSRMatrixData( MatrixData &A )
{
    auto ptr = dynamic_cast<CSRMatrixData<Policy> *>( &A );
    AMP_INSIST( ptr, "dynamic cast from const MatrixData to const CSRMatrixData failed" );
    return ptr;
}


template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::mult( std::shared_ptr<const Vector> in,
                                               MatrixData const &A,
                                               std::shared_ptr<Vector> out )
{
    AMP_ASSERT( in && out );

    using gidx_t   = typename Policy::gidx_t;
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );

    auto [nnz, cols, coeffs] = csrData->getCSRData();

    auto memType = AMP::Utilities::getMemoryType( cols );
    AMP_INSIST( memType == AMP::Utilities::MemoryType::host ||
                    memType == AMP::Utilities::MemoryType::unregistered,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    const auto nRows = csrData->numLocalRows();
    auto maxColLen   = *std::max_element( nnz, nnz + nRows );

    std::vector<size_t> rcols( maxColLen );
    std::vector<scalar_t> vvals( maxColLen );

    auto beginRow = csrData->beginRow();

    lidx_t offset = 0;
    for ( lidx_t row = 0; row < nRows; ++row ) {

        const auto nCols = nnz[row];

        const auto cloc = &cols[offset];
        const auto vloc = &coeffs[offset];

        std::transform(
            cloc, cloc + nCols, rcols.begin(), []( gidx_t col ) -> size_t { return col; } );

        in->getValuesByGlobalID( nCols, rcols.data(), vvals.data() );

        scalar_t val =
            std::inner_product( vloc, vloc + nCols, vvals.data(), static_cast<scalar_t>( 0.0 ) );

        out->setValueByGlobalID( static_cast<size_t>( beginRow + row ), val );

        offset += nCols;
    }

    out->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::multTranspose( std::shared_ptr<const Vector> in,
                                                        MatrixData const &A,
                                                        std::shared_ptr<Vector> out )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::scale( AMP::Scalar alpha_in, MatrixData &A )
{
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );

    auto [nnz_per_row, cols, coeffs] = csrData->getCSRData();

    auto memType = AMP::Utilities::getMemoryType( cols );
    AMP_INSIST( memType == AMP::Utilities::MemoryType::host ||
                    memType == AMP::Utilities::MemoryType::unregistered,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    const auto nRows = csrData->numLocalRows();
    const auto nnz = std::accumulate( nnz_per_row, nnz_per_row + nRows, static_cast<lidx_t>( 0 ) );

    auto alpha = static_cast<scalar_t>( alpha_in );
    std::transform( const_cast<scalar_t *>( coeffs ),
                    const_cast<scalar_t *>( coeffs ) + nnz,
                    const_cast<scalar_t *>( coeffs ),
                    [=]( scalar_t val ) { return alpha * val; } );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::matMultiply( MatrixData const &Am,
                                                      MatrixData const &Bm,
                                                      MatrixData &Cm )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::axpy( AMP::Scalar alpha_in,
                                               const MatrixData &X,
                                               MatrixData &Y )
{
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    const auto csrDataX = getCSRMatrixData<Policy>( const_cast<MatrixData &>( X ) );
    const auto [nnz_x, cols_x, coeffs_x] = csrDataX->getCSRData();

    auto csrDataY                  = getCSRMatrixData<Policy>( Y );
    auto [nnz_y, cols_y, coeffs_y] = csrDataY->getCSRData();

    auto memType_x = AMP::Utilities::getMemoryType( cols_x );
    AMP_INSIST( memType_x == AMP::Utilities::MemoryType::host ||
                    memType_x == AMP::Utilities::MemoryType::unregistered,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    auto memType_y = AMP::Utilities::getMemoryType( cols_y );
    AMP_INSIST( memType_y == AMP::Utilities::MemoryType::host ||
                    memType_y == AMP::Utilities::MemoryType::unregistered,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    const auto nRows = csrDataX->numLocalRows();
    const auto nnz   = std::accumulate( nnz_x, nnz_x + nRows, static_cast<lidx_t>( 0 ) );

    auto alpha = static_cast<scalar_t>( alpha_in );
    auto y_p   = const_cast<scalar_t *>( coeffs_y );
    for ( lidx_t i = 0; i < nnz; ++i )
        y_p[i] += alpha * coeffs_x[i];
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::setScalar( AMP::Scalar alpha_in, MatrixData &A )
{
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );

    auto [nnz_per_row, cols, coeffs] = csrData->getCSRData();

    auto memType = AMP::Utilities::getMemoryType( cols );
    AMP_INSIST( memType == AMP::Utilities::MemoryType::host ||
                    memType == AMP::Utilities::MemoryType::unregistered,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    const auto nRows = csrData->numLocalRows();
    const auto nnz = std::accumulate( nnz_per_row, nnz_per_row + nRows, static_cast<lidx_t>( 0 ) );

    auto alpha = static_cast<scalar_t>( alpha_in );
    std::fill( const_cast<scalar_t *>( coeffs ), const_cast<scalar_t *>( coeffs ) + nnz, alpha );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::zero( MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;
    setScalar( static_cast<scalar_t>( 0.0 ), A );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::setDiagonal( std::shared_ptr<const Vector> in,
                                                      MatrixData &A )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    // constrain to one data block for now
    AMP_ASSERT( in && in->numberOfDataBlocks() == 1 && in->isType<scalar_t>( 0 ) );

    const scalar_t *vvals_p = in->getRawDataBlock<scalar_t>();

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );

    auto [nnz, cols, coeffs] = csrData->getCSRData();

    auto memType = AMP::Utilities::getMemoryType( cols );
    AMP_INSIST( memType == AMP::Utilities::MemoryType::host ||
                    memType == AMP::Utilities::MemoryType::unregistered,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    const auto nRows = csrData->numLocalRows();

    auto beginRow = csrData->beginRow();

    auto vals_p = const_cast<scalar_t *>( coeffs );

    lidx_t offset = 0;
    for ( lidx_t row = 0; row < nRows; ++row ) {
        const auto ncols = nnz[row];
        for ( lidx_t icol = 0; icol < ncols; ++icol ) {
            if ( cols[offset + icol] == static_cast<gidx_t>( beginRow + row ) ) {
                vals_p[offset + icol] = vvals_p[row];
                break;
            }
        }
        offset += nnz[row];
    }
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::setIdentity( MatrixData &A )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    zero( A );

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );

    auto [nnz, cols, coeffs] = csrData->getCSRData();

    auto memType = AMP::Utilities::getMemoryType( cols );
    AMP_INSIST( memType == AMP::Utilities::MemoryType::host ||
                    memType == AMP::Utilities::MemoryType::unregistered,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    const auto nRows = csrData->numLocalRows();

    auto beginRow = csrData->beginRow();

    auto vals_p = const_cast<scalar_t *>( coeffs );

    lidx_t offset = 0;
    for ( lidx_t row = 0; row < nRows; ++row ) {
        const auto ncols = nnz[row];
        for ( lidx_t icol = 0; icol < ncols; ++icol ) {
            if ( cols[offset + icol] == static_cast<gidx_t>( beginRow + row ) ) {
                vals_p[offset + icol] = static_cast<scalar_t>( 1.0 );
                break;
            }
        }
        offset += nnz[row];
    }
}

template<typename Policy>
AMP::Scalar CSRMatrixOperationsDefault<Policy>::L1Norm( MatrixData const &A ) const
{
    AMP_ERROR( "Not implemented" );
}

} // namespace AMP::LinearAlgebra
