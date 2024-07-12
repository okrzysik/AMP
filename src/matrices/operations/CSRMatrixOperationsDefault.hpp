#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/CSRMatrixOperationsDefault.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>

#include "ProfilerApp.h"

namespace AMP::LinearAlgebra {

template<typename Policy, typename Allocator>
void CSRMatrixOperationsDefault<Policy,Allocator>::mult( std::shared_ptr<const Vector> in,
                                               MatrixData const &A,
                                               std::shared_ptr<Vector> out )
{
    PROFILE( "CSRMatrixOperationsDefault::mult" );
    AMP_DEBUG_ASSERT( in && out );
    AMP_DEBUG_ASSERT( in->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::UNCHANGED );

    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy,Allocator>( const_cast<MatrixData &>( A ) );

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

template<typename Policy,typename Allocator>
void CSRMatrixOperationsDefault<Policy,Allocator>::multTranspose( std::shared_ptr<const Vector> in,
                                                        MatrixData const &A,
                                                        std::shared_ptr<Vector> out )
{
    PROFILE( "CSRMatrixOperationsDefault::multTranspose" );
    
    // this is not meant to be an optimized version. It is provided for completeness
    AMP_DEBUG_ASSERT( in && out );

    out->zero();

    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy,Allocator>( const_cast<MatrixData &>( A ) );
    
    AMP_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    auto inData                 = in->getVectorData();
    const scalar_t *inDataBlock = inData->getRawDataBlock<scalar_t>( 0 );

    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );

    {
        PROFILE( "CSRMatrixOperationsDefault::multTranspose (d)" );
	auto [nnz, cols, cols_loc, coeffs] = csrData->getCSRDiagData();
	const auto num_unq = csrData->numLocalColumns();
	
        std::vector<scalar_t> vvals( num_unq, 0.0 );
	std::vector<size_t> rcols( num_unq );
	
        lidx_t offset = 0;
        for ( lidx_t row = 0; row < nRows; ++row ) {

            const auto ncols = nnz[row];
            const auto cloc = &cols_loc[offset];
            const auto vloc = &coeffs[offset];
            const auto val = inDataBlock[ row ];

            for ( lidx_t j = 0; j < ncols; ++j ) {
	        rcols[cloc[j]] = cols[offset + j];
                vvals[cloc[j]] += vloc[j] * val;
            }

            offset += ncols;
        }

	// Write out data, adding to any already present
	out->addValuesByGlobalID( num_unq, rcols.data(), vvals.data() );
    }
    out->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );

    if ( csrData->hasOffDiag() ) {
        PROFILE( "CSRMatrixOperationsDefault::multTranspose (od)" );
	auto [nnz, cols, cols_loc, coeffs] = csrData->getCSROffDiagData();
	
	std::vector<size_t> rcols;
	csrData->getOffDiagColumnMap( rcols );
	const auto num_unq = rcols.size();
	
        std::vector<scalar_t> vvals( num_unq, 0.0 );
	
        lidx_t offset = 0;
        for ( lidx_t row = 0; row < nRows; ++row ) {

            const auto ncols = nnz[row];
	    if ( ncols == 0 ) { continue; }
	    
            const auto cloc = &cols_loc[offset];
            const auto vloc = &coeffs[offset];
            const auto val = inDataBlock[ row ];

            for ( lidx_t j = 0; j < ncols; ++j ) {
                vvals[cloc[j]] += vloc[j] * val;
            }

            offset += ncols;
        }

	// convert colmap to size_t and write out data
	out->addValuesByGlobalID( num_unq, rcols.data(), vvals.data() );
    }
}

template<typename Policy, typename Allocator>
void CSRMatrixOperationsDefault<Policy,Allocator>::scale( AMP::Scalar alpha_in, MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy,Allocator>( const_cast<MatrixData &>( A ) );

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

template<typename Policy, typename Allocator>
void CSRMatrixOperationsDefault<Policy,Allocator>::matMultiply( MatrixData const &,
                                                      MatrixData const &,
                                                      MatrixData & )
{
    AMP_WARNING( "SpGEMM for CSRMatrixOperationsDefault not implemented" );
}

template<typename Policy, typename Allocator>
void CSRMatrixOperationsDefault<Policy,Allocator>::axpy( AMP::Scalar alpha_in,
                                               const MatrixData &X,
                                               MatrixData &Y )
{
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    const auto csrDataX = getCSRMatrixData<Policy,Allocator>( const_cast<MatrixData &>( X ) );
    const auto [nnz_d_x, cols_d_x, cols_loc_d_x, coeffs_d_x] = csrDataX->getCSRDiagData();

    auto csrDataY                                      = getCSRMatrixData<Policy,Allocator>( Y );
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
        for ( lidx_t i = 0; i < tnnz; ++i ) {
            y_p[i] += alpha * coeffs_d_x[i];
        }
    }

    if ( csrDataX->hasOffDiag() ) {
        const auto [nnz_od_x, cols_od_x, cols_loc_od_x, coeffs_od_x] =
            csrDataX->getCSROffDiagData();
        auto [nnz_od_y, cols_od_y, cols_loc_od_y, coeffs_od_y] = csrDataY->getCSROffDiagData();
        const auto tnnz = csrDataX->numberOfNonZerosOffDiag();
        auto y_p        = const_cast<scalar_t *>( coeffs_od_y );
        for ( lidx_t i = 0; i < tnnz; ++i ) {
            y_p[i] += alpha * coeffs_od_x[i];
        }
    }
}

template<typename Policy, typename Allocator>
void CSRMatrixOperationsDefault<Policy,Allocator>::setScalar( AMP::Scalar alpha_in, MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy,Allocator>( const_cast<MatrixData &>( A ) );

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

template<typename Policy, typename Allocator>
void CSRMatrixOperationsDefault<Policy,Allocator>::zero( MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;
    setScalar( static_cast<scalar_t>( 0.0 ), A );
}

template<typename Policy, typename Allocator>
void CSRMatrixOperationsDefault<Policy,Allocator>::setDiagonal( std::shared_ptr<const Vector> in,
                                                      MatrixData &A )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    // constrain to one data block for now
    AMP_DEBUG_ASSERT( in && in->numberOfDataBlocks() == 1 && in->isType<scalar_t>( 0 ) );

    const scalar_t *vvals_p = in->getRawDataBlock<scalar_t>();

    auto csrData = getCSRMatrixData<Policy,Allocator>( const_cast<MatrixData &>( A ) );

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

template<typename Policy, typename Allocator>
void CSRMatrixOperationsDefault<Policy,Allocator>::setIdentity( MatrixData &A )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    zero( A );

    auto csrData = getCSRMatrixData<Policy,Allocator>( const_cast<MatrixData &>( A ) );

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

template<typename Policy, typename Allocator>
AMP::Scalar CSRMatrixOperationsDefault<Policy,Allocator>::L1Norm( MatrixData const &A ) const
{
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy,Allocator>( const_cast<MatrixData &>( A ) );

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
