#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/CSRMatrixOperationsDefault.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Vector.h"

#include <algorithm>

#include "ProfilerApp.h"

namespace AMP::LinearAlgebra {

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixOperationsDefault<Policy, Allocator>::mult( std::shared_ptr<const Vector> in,
                                                          MatrixData const &A,
                                                          std::shared_ptr<Vector> out )
{
    PROFILE( "CSRMatrixOperationsDefault::mult" );
    AMP_DEBUG_ASSERT( in && out );
    AMP_DEBUG_ASSERT( in->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::UNCHANGED );

    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>(
        const_cast<MatrixData &>( A ) );

    // this ensures that the separate on/off diag calls operator the same way
    // e.g. by accumulating into output data
    out->zero();

    auto inData                 = in->getVectorData();
    const scalar_t *inDataBlock = inData->getRawDataBlock<scalar_t>( 0 );
    const auto &ghosts          = inData->getGhosts();
    auto outData                = out->getVectorData();
    scalar_t *outDataBlock      = outData->getRawDataBlock<scalar_t>( 0 );

    AMP_INSIST( csrData->d_memory_location != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is not implemented for device memory" );

    AMP_INSIST(
        1 == inData->numberOfDataBlocks(),
        "CSRMatrixOperationsDefault::mult only implemented for vectors with one data block" );

    AMP_INSIST(
        ghosts.size() == inData->getGhostSize(),
        "CSRMatrixOperationsDefault::mult only implemented for vectors with accessible ghosts" );

    AMP_ASSERT( inDataBlock && outDataBlock );

    {
        PROFILE( "CSRMatrixOperationsDefault::mult (local)" );
        d_localops_diag.mult( inDataBlock, csrData->getDiagMatrix(), outDataBlock );
    }

    if ( csrData->hasOffDiag() ) {
        PROFILE( "CSRMatrixOperationsDefault::mult (ghost)" );
        d_localops_offd.mult( inDataBlock, csrData->getOffdMatrix(), outDataBlock );
    }
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixOperationsDefault<Policy, Allocator>::multTranspose( std::shared_ptr<const Vector> in,
                                                                   MatrixData const &A,
                                                                   std::shared_ptr<Vector> out )
{
    PROFILE( "CSRMatrixOperationsDefault::multTranspose" );

    // this is not meant to be an optimized version. It is provided for completeness
    AMP_DEBUG_ASSERT( in && out );

    out->zero();

    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>(
        const_cast<MatrixData &>( A ) );

    AMP_INSIST( csrData->d_memory_location != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is not implemented for device memory" );

    auto inData                 = in->getVectorData();
    const scalar_t *inDataBlock = inData->getRawDataBlock<scalar_t>( 0 );

    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );

    {
        PROFILE( "CSRMatrixOperationsDefault::multTranspose (d)" );

        std::vector<scalar_t> vvals( num_unq, 0.0 );
        std::vector<size_t> rcols( num_unq );

        d_localops_diag.multTranspose( inDataBlock, csrData->getDiagMatrix(), vvals, rcols );

        // Write out data, adding to any already present
        out->addValuesByGlobalID( num_unq, rcols.data(), vvals.data() );
    }
    out->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );

    if ( csrData->hasOffDiag() ) {
        PROFILE( "CSRMatrixOperationsDefault::multTranspose (od)" );
        d_localops_offd.multTranspose( inDataBlock, csrData->getOffdMatrix(), SOMETHING );
        // write out data
        out->addValuesByGlobalID( num_unq, rcols.data(), vvals.data() );
    }
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixOperationsDefault<Policy, Allocator>::scale( AMP::Scalar alpha_in, MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>(
        const_cast<MatrixData &>( A ) );

    AMP_INSIST( csrData->d_memory_location != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is not implemented for device memory" );

    auto alpha = static_cast<scalar_t>( alpha_in );
    d_localops_diag.scale( alpha, csrData->getDiagMatrix() );
    if ( csrData->hasOffDiag() ) {
        d_localops_offd.scale( alpha, csrData->getOffdMatrix() );
    }
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixOperationsDefault<Policy, Allocator>::matMultiply( MatrixData const &,
                                                                 MatrixData const &,
                                                                 MatrixData & )
{
    AMP_WARNING( "SpGEMM for CSRMatrixOperationsDefault not implemented" );
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixOperationsDefault<Policy, Allocator>::axpy( AMP::Scalar alpha_in,
                                                          const MatrixData &X,
                                                          MatrixData &Y )
{
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    const auto csrDataX = getCSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>(
        const_cast<MatrixData &>( X ) );
    auto csrDataY = getCSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>( Y );

    AMP_INSIST( csrDataX->d_memory_location != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is not implemented for device memory" );

    AMP_INSIST( csrDataY->d_memory_location != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is not implemented for device memory" );

    auto alpha = static_cast<scalar_t>( alpha_in );
    d_localops_diag.axpy( alpha, csrDataX->getDiagMatrix(), csrDataY->getDiagMatrix() );
    if ( csrDataX->hasOffDiag() ) {
        d_localops_offd.axpy( alpha, csrDataX->getOffdMatrix(), csrDataY->getOffdMatrix() );
    }
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixOperationsDefault<Policy, Allocator>::setScalar( AMP::Scalar alpha_in, MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>(
        const_cast<MatrixData &>( A ) );

    AMP_INSIST( csrData->d_memory_location != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is not implemented for device memory" );

    auto alpha = static_cast<scalar_t>( alpha_in );
    d_localops_diag.setScalar( alpha, csrData->getDiagMatrix() );
    if ( csrData->hasOffDiag() ) {
        d_localops_offd.setScalar( alpha, csrData->getOffdMatrix() );
    }
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixOperationsDefault<Policy, Allocator>::zero( MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;
    setScalar( static_cast<scalar_t>( 0.0 ), A );
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixOperationsDefault<Policy, Allocator>::setDiagonal( std::shared_ptr<const Vector> in,
                                                                 MatrixData &A )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    // constrain to one data block for now
    AMP_DEBUG_ASSERT( in && in->numberOfDataBlocks() == 1 && in->isType<scalar_t>( 0 ) );

    const scalar_t *vvals_p = in->getRawDataBlock<scalar_t>();

    auto csrData = getCSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>(
        const_cast<MatrixData &>( A ) );

    AMP_INSIST( csrData->d_memory_location != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is not implemented for device memory" );

    d_localops_diag.setDiagonal( vvals_p, csrData->getDiagMatrix() );
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixOperationsDefault<Policy, Allocator>::setIdentity( MatrixData &A )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>(
        const_cast<MatrixData &>( A ) );

    AMP_INSIST( csrData->d_memory_location != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is not implemented for device memory" );

    zero( A );

    d_localops_diag.setIdentity( csrData->getDiagMatrix() );
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixOperationsDefault<Policy, Allocator>::extractDiagonal( MatrixData const &A,
                                                                     std::shared_ptr<Vector> buf )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>(
        const_cast<MatrixData &>( A ) );

    AMP_ASSERT( buf && buf->numberOfDataBlocks() == 1 );
    AMP_ASSERT( buf->isType<scalar_t>( 0 ) );

    auto *rawVecData = buf->getRawDataBlock<scalar_t>();
    auto memTypeV    = AMP::Utilities::getMemoryType( rawVecData );
    AMP_INSIST( memTypeV < AMP::Utilities::MemoryType::device &&
                    csrData->d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault::extractDiagonal not implemented for device memory" );

    d_localops_diag.extractDiagonal( csrData->getDiagMatrix(), rawVecData );
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
AMP::Scalar CSRMatrixOperationsDefault<Policy, Allocator>::LinfNorm( MatrixData const &A ) const
{
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>(
        const_cast<MatrixData &>( A ) );

    AMP_INSIST( csrData->d_memory_location != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is not implemented for device memory" );

    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );
    std::vector<scalar_t> rowSums( nRows, 0.0 );

    d_localops_diag.LinfNorm( csrData->getDiagMatrix(), rowSums.data() );
    if ( csrData->hasOffDiag() ) {
        d_localops_offd.LinfNorm( csrData->getOffdMatrix(), rowSums.data() );
    }

    // Reduce row sums to get global Linf norm
    auto max_norm = *std::max_element( rowSums.begin(), rowSums.end() );
    AMP_MPI comm  = csrData->getComm();
    return comm.maxReduce<scalar_t>( max_norm );
}

} // namespace AMP::LinearAlgebra
