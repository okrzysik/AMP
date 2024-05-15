#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/CSRMatrixOperationsDefault.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>

#include "ProfilerApp.h"

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
    PROFILE( "CSRMatrixOperationsDefault::mult" );
    AMP_DEBUG_ASSERT( in && out );
    AMP_DEBUG_ASSERT( in->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::UNCHANGED );

    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );

    auto [nnz_d, cols_d, cols_loc_d, coeffs_d] = csrData->getCSRDiagData();

    auto inData                 = in->getVectorData();
    const scalar_t *inDataBlock = inData->getRawDataBlock<scalar_t>( 0 );
    const auto &ghosts          = inData->getGhosts();
    auto outData                = out->getVectorData();
    scalar_t *outDataBlock      = outData->getRawDataBlock<scalar_t>( 0 );

    AMP_INSIST( AMP::Utilities::getMemoryType( cols_loc_d ) != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    AMP_INSIST(
        1 == inData->numberOfDataBlocks(),
        "CSRMatrixOperationsDefault::mult only implemented for vectors with one data block" );

    AMP_INSIST(
        ghosts.size() == inData->getGhostSize(),
        "CSRMatrixOperationsDefault::mult only implemented for vectors with accessible ghosts" );

    AMP_ASSERT( inDataBlock && outDataBlock );

    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );

    {
        PROFILE( "CSRMatrixOperationsDefault::mult (local)" );
        lidx_t offset = 0;
        for ( lidx_t row = 0; row < nRows; ++row ) {
            const auto nCols = nnz_d[row];
            const auto cloc  = &cols_loc_d[offset];
            const auto vloc  = &coeffs_d[offset];

            outDataBlock[row] = 0.0;
            for ( lidx_t c = 0; c < nCols; ++c ) {
                outDataBlock[row] += vloc[c] * inDataBlock[cloc[c]];
            }

            offset += nCols;
        }
    }

    if ( csrData->hasOffDiag() ) {
        PROFILE( "CSRMatrixOperationsDefault::mult (ghost)" );
        auto [nnz_od, cols_od, cols_loc_od, coeffs_od] = csrData->getCSROffDiagData();
        lidx_t offset                                  = 0;

        for ( lidx_t row = 0; row < nRows; ++row ) {
            const auto nCols = nnz_od[row];
            const auto cloc  = &cols_loc_od[offset];
            const auto vloc  = &coeffs_od[offset];

            for ( lidx_t c = 0; c < nCols; ++c ) {
                outDataBlock[row] += vloc[c] * ghosts[cloc[c]];
            }

            offset += nCols;
        }
    }
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::multTranspose( std::shared_ptr<const Vector> in,
                                                        MatrixData const &A,
                                                        std::shared_ptr<Vector> out )
{
    // this is not meant to be an optimized version. It is provided for completeness
    AMP_DEBUG_ASSERT( in && out );
    AMP_DEBUG_ASSERT( in->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::UNCHANGED );

    out->zero();

    using gidx_t   = typename Policy::gidx_t;
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );

    auto [nnz_d, cols_d, cols_loc_d, coeffs_d] = csrData->getCSRDiagData();

    auto memType = AMP::Utilities::getMemoryType( cols_loc_d );
    AMP_INSIST( memType != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );
    auto maxColLen   = *std::max_element( nnz_d, nnz_d + nRows );
    std::vector<size_t> rcols( maxColLen );
    std::vector<scalar_t> vvals( maxColLen );

    {
        auto maxColLen = *std::max_element( nnz_d, nnz_d + nRows );
        std::vector<size_t> rcols( maxColLen );
        std::vector<scalar_t> vvals( maxColLen );
        lidx_t offset = 0;
        for ( lidx_t row = 0; row < nRows; ++row ) {

            const auto nCols = nnz_d[row];

            const auto cloc = &cols_d[offset];
            const auto vloc = &coeffs_d[offset];

            std::transform(
                cloc, cloc + nCols, rcols.begin(), []( gidx_t col ) -> size_t { return col; } );

            const auto val = in->getValueByGlobalID( row );

            for ( lidx_t icol = 0; icol < nCols; ++icol ) {
                vvals[icol] = vloc[icol] * val;
            }

            out->addValuesByGlobalID( nCols, rcols.data(), vvals.data() );

            offset += nCols;
        }
    }

    if ( csrData->hasOffDiag() ) {
        auto [nnz_od, cols_od, cols_loc_od, coeffs_od] = csrData->getCSRDiagData();
        auto maxColLen = *std::max_element( nnz_od, nnz_od + nRows );
        std::vector<size_t> rcols( maxColLen );
        std::vector<scalar_t> vvals( maxColLen );
        lidx_t offset = 0;
        for ( lidx_t row = 0; row < nRows; ++row ) {

            const auto nCols = nnz_od[row];

            const auto cloc = &cols_od[offset];
            const auto vloc = &coeffs_od[offset];

            std::transform(
                cloc, cloc + nCols, rcols.begin(), []( gidx_t col ) -> size_t { return col; } );

            const auto val = in->getValueByGlobalID( row );

            for ( lidx_t icol = 0; icol < nCols; ++icol ) {
                vvals[icol] = vloc[icol] * val;
            }

            out->addValuesByGlobalID( nCols, rcols.data(), vvals.data() );

            offset += nCols;
        }
    }
    // consistent add because some values might be remote
    out->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::scale( AMP::Scalar alpha_in, MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );

    auto [nnz_d, cols_d, cols_loc_d, coeffs_d] = csrData->getCSRDiagData();

    auto memType = AMP::Utilities::getMemoryType( cols_loc_d );
    AMP_INSIST( memType != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    const auto tnnz_d = csrData->numberOfNonZerosDiag();

    auto alpha = static_cast<scalar_t>( alpha_in );
    std::transform( const_cast<scalar_t *>( coeffs_d ),
                    const_cast<scalar_t *>( coeffs_d ) + tnnz_d,
                    const_cast<scalar_t *>( coeffs_d ),
                    [=]( scalar_t val ) { return alpha * val; } );

    if ( csrData->hasOffDiag() ) {
        const auto tnnz_od                             = csrData->numberOfNonZerosOffDiag();
        auto [nnz_od, cols_od, cols_loc_od, coeffs_od] = csrData->getCSROffDiagData();
        std::transform( const_cast<scalar_t *>( coeffs_od ),
                        const_cast<scalar_t *>( coeffs_od ) + tnnz_od,
                        const_cast<scalar_t *>( coeffs_od ),
                        [=]( scalar_t val ) { return alpha * val; } );
    }
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::matMultiply( MatrixData const &,
                                                      MatrixData const &,
                                                      MatrixData & )
{
    AMP_WARNING( "SpGEMM for CSRMatrixOperationsDefault not implemented" );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::axpy( AMP::Scalar alpha_in,
                                               const MatrixData &X,
                                               MatrixData &Y )
{
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    const auto csrDataX = getCSRMatrixData<Policy>( const_cast<MatrixData &>( X ) );
    const auto [nnz_d_x, cols_d_x, cols_loc_d_x, coeffs_d_x] = csrDataX->getCSRDiagData();

    auto csrDataY                                      = getCSRMatrixData<Policy>( Y );
    auto [nnz_d_y, cols_d_y, cols_loc_d_y, coeffs_d_y] = csrDataY->getCSRDiagData();

    auto memType_x = AMP::Utilities::getMemoryType( cols_loc_d_x );
    AMP_INSIST( memType_x != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    auto memType_y = AMP::Utilities::getMemoryType( cols_loc_d_y );
    AMP_INSIST( memType_y != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    auto alpha = static_cast<scalar_t>( alpha_in );

    {
        const auto tnnz = csrDataX->numberOfNonZerosDiag();
        auto y_p        = const_cast<scalar_t *>( coeffs_d_y );
        for ( gidx_t i = 0; i < tnnz; ++i ) {
            y_p[i] += alpha * coeffs_d_x[i];
        }
    }

    if ( csrDataX->hasOffDiag() ) {
        const auto [nnz_od_x, cols_od_x, cols_loc_od_x, coeffs_od_x] =
            csrDataX->getCSROffDiagData();
        auto [nnz_od_y, cols_od_y, cols_loc_od_y, coeffs_od_y] = csrDataY->getCSROffDiagData();
        const auto tnnz = csrDataX->numberOfNonZerosOffDiag();
        auto y_p        = const_cast<scalar_t *>( coeffs_od_y );
        for ( gidx_t i = 0; i < tnnz; ++i ) {
            y_p[i] += alpha * coeffs_od_x[i];
        }
    }
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::setScalar( AMP::Scalar alpha_in, MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );

    auto [nnz_d, cols_d, cols_loc_d, coeffs_d] = csrData->getCSRDiagData();

    auto memType = AMP::Utilities::getMemoryType( cols_loc_d );
    AMP_INSIST( memType != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    const auto tnnz_d = csrData->numberOfNonZerosDiag();

    auto alpha = static_cast<scalar_t>( alpha_in );
    std::fill(
        const_cast<scalar_t *>( coeffs_d ), const_cast<scalar_t *>( coeffs_d ) + tnnz_d, alpha );

    if ( csrData->hasOffDiag() ) {
        const auto tnnz_od                             = csrData->numberOfNonZerosOffDiag();
        auto [nnz_od, cols_od, cols_loc_od, coeffs_od] = csrData->getCSROffDiagData();
        std::fill( const_cast<scalar_t *>( coeffs_od ),
                   const_cast<scalar_t *>( coeffs_od ) + tnnz_od,
                   alpha );
    }
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
    AMP_DEBUG_ASSERT( in && in->numberOfDataBlocks() == 1 && in->isType<scalar_t>( 0 ) );

    const scalar_t *vvals_p = in->getRawDataBlock<scalar_t>();

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );

    auto [nnz_d, cols_d, cols_loc_d, coeffs_d] = csrData->getCSRDiagData();

    auto memType = AMP::Utilities::getMemoryType( cols_loc_d );
    AMP_INSIST( memType != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );

    auto beginRow = csrData->beginRow();

    auto vals_p = const_cast<scalar_t *>( coeffs_d );

    lidx_t offset = 0;
    for ( lidx_t row = 0; row < nRows; ++row ) {
        const auto ncols = nnz_d[row];
        for ( lidx_t icol = 0; icol < ncols; ++icol ) {
            if ( cols_d[offset + icol] == static_cast<gidx_t>( beginRow + row ) ) {
                vals_p[offset + icol] = vvals_p[row];
                break;
            }
        }
        offset += nnz_d[row];
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

    auto [nnz_d, cols_d, cols_loc_d, coeffs_d] = csrData->getCSRDiagData();

    auto memType = AMP::Utilities::getMemoryType( cols_loc_d );
    AMP_INSIST( memType != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );

    auto beginRow = csrData->beginRow();

    auto vals_p = const_cast<scalar_t *>( coeffs_d );

    lidx_t offset = 0;
    for ( lidx_t row = 0; row < nRows; ++row ) {
        const auto ncols = nnz_d[row];
        for ( lidx_t icol = 0; icol < ncols; ++icol ) {
            if ( cols_d[offset + icol] == static_cast<gidx_t>( beginRow + row ) ) {
                vals_p[offset + icol] = static_cast<scalar_t>( 1.0 );
                break;
            }
        }
        offset += nnz_d[row];
    }
}

template<typename Policy>
AMP::Scalar CSRMatrixOperationsDefault<Policy>::L1Norm( MatrixData const &A ) const
{
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );

    auto [nnz_d, cols_d, cols_loc_d, coeffs_d] = csrData->getCSRDiagData();

    auto memType = AMP::Utilities::getMemoryType( cols_loc_d );
    AMP_INSIST( memType != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    const auto ncols = csrData->numGlobalColumns();
    std::vector<scalar_t> col_norms( ncols, 0.0 );

    const size_t tnnz_d = static_cast<size_t>( csrData->numberOfNonZerosDiag() );
    for ( size_t i = 0; i < tnnz_d; ++i ) {
        col_norms[cols_d[i]] += std::abs( coeffs_d[i] );
    }

    if ( csrData->hasOffDiag() ) {
        const size_t tnnz_od = static_cast<size_t>( csrData->numberOfNonZerosOffDiag() );
        auto [nnz_od, cols_od, cols_loc_od, coeffs_od] = csrData->getCSROffDiagData();
        for ( size_t i = 0; i < tnnz_od; ++i ) {
            col_norms[cols_od[i]] += std::abs( coeffs_od[i] );
        }
    }
    // Reduce partial column sums across all ranks to get full column norms
    AMP_MPI comm = csrData->getComm();
    comm.sumReduce<scalar_t>( col_norms.data(), ncols );

    return *std::max_element( col_norms.begin(), col_norms.end() );
}

} // namespace AMP::LinearAlgebra
