#ifndef included_CSRMatrixOperationsDevice_HPP_
#define included_CSRMatrixOperationsDevice_HPP_

#include "AMP/matrices/CSRPolicy.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/device/CSRLocalMatrixOperationsDevice.h"
#include "AMP/matrices/operations/device/CSRMatrixOperationsDevice.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/typeid.h"

#include "AMP/utils/memory.h"

#include "thrust/device_vector.h"
#include "thrust/execution_policy.h"
#include "thrust/extrema.h"

#include <algorithm>

#include "ProfilerApp.h"

namespace AMP::LinearAlgebra {

template<typename Policy, class Allocator>
void CSRMatrixOperationsDevice<Policy, Allocator>::mult( std::shared_ptr<const Vector> in,
                                                         MatrixData const &A,
                                                         std::shared_ptr<Vector> out )
{
    PROFILE( "CSRMatrixOperationsDevice::mult" );
    AMP_DEBUG_ASSERT( in && out );
    AMP_DEBUG_ASSERT( in->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::UNCHANGED );

    auto csrData = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );

    AMP_DEBUG_ASSERT( csrData );

    auto diagMatrix = csrData->getDiagMatrix();
    auto offdMatrix = csrData->getOffdMatrix();

    AMP_DEBUG_ASSERT( diagMatrix && offdMatrix );

    out->zero();

    auto inData                 = in->getVectorData();
    const scalar_t *inDataBlock = inData->getRawDataBlock<scalar_t>( 0 );
    auto outData                = out->getVectorData();
    scalar_t *outDataBlock      = outData->getRawDataBlock<scalar_t>( 0 );

    AMP_DEBUG_INSIST( csrData->d_memory_location != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsDevice is not implemented for device memory" );

    AMP_DEBUG_INSIST( csrData->d_memory_location == AMP::Utilities::getMemoryType( inDataBlock ),
                      "Input vector from wrong memory space" );

    AMP_DEBUG_INSIST( csrData->d_memory_location == AMP::Utilities::getMemoryType( outDataBlock ),
                      "Output vector from wrong memory space" );

    AMP_DEBUG_INSIST(
        1 == inData->numberOfDataBlocks(),
        "CSRMatrixOperationsDevice::mult only implemented for vectors with one data block" );

    AMP_ASSERT( inDataBlock && outDataBlock );

    {
        PROFILE( "CSRMatrixOperationsDevice::mult(local)" );
        CSRLocalMatrixOperationsDevice<Policy, Allocator>::mult(
            inDataBlock, diagMatrix, outDataBlock );
    }

    if ( csrData->hasOffDiag() ) {
        PROFILE( "CSRMatrixOperationsDevice::mult(ghost)" );
        using scalarAllocator_t =
            typename std::allocator_traits<Allocator>::template rebind_alloc<scalar_t>;
        const auto nGhosts = offdMatrix->numUniqueColumns();
        scalarAllocator_t alloc;
        scalar_t *ghosts = alloc.allocate( nGhosts );
        if constexpr ( std::is_same_v<size_t, gidx_t> ) {
            // column map can be passed to get ghosts function directly
            size_t *colMap = offdMatrix->getColumnMap();
            in->getGhostValuesByGlobalID( nGhosts, colMap, ghosts );
        } else {
            std::vector<size_t> colMap;
            offdMatrix->getColumnMap( colMap );
            in->getGhostValuesByGlobalID( nGhosts, colMap.data(), ghosts );
        }

        CSRLocalMatrixOperationsDevice<Policy, Allocator>::mult( ghosts, offdMatrix, outDataBlock );
        alloc.deallocate( ghosts, nGhosts );
    }
}

template<typename Policy, class Allocator>
void CSRMatrixOperationsDevice<Policy, Allocator>::multTranspose( std::shared_ptr<const Vector> in,
                                                                  MatrixData const &A,
                                                                  std::shared_ptr<Vector> out )
{
    AMP_WARNING( "multTranspose not enabled for device." );
}

template<typename Policy, class Allocator>
void CSRMatrixOperationsDevice<Policy, Allocator>::scale( AMP::Scalar alpha_in, MatrixData &A )
{
    auto csrData = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );

    AMP_DEBUG_ASSERT( csrData );

    auto diagMatrix = csrData->getDiagMatrix();
    AMP_DEBUG_ASSERT( diagMatrix );

    auto alpha = static_cast<scalar_t>( alpha_in );
    CSRLocalMatrixOperationsDevice<Policy, Allocator>::scale( alpha, diagMatrix );

    if ( csrData->hasOffDiag() ) {
        auto offdMatrix = csrData->getOffdMatrix();
        AMP_DEBUG_ASSERT( offdMatrix );
        CSRLocalMatrixOperationsDevice<Policy, Allocator>::scale( alpha, offdMatrix );
    }
}

template<typename Policy, class Allocator>
void CSRMatrixOperationsDevice<Policy, Allocator>::matMatMult( std::shared_ptr<MatrixData>,
                                                               std::shared_ptr<MatrixData>,
                                                               std::shared_ptr<MatrixData> )
{
    AMP_WARNING( "matMatMult for CSRMatrixOperationsDevice not implemented" );
}

template<typename Policy, class Allocator>
void CSRMatrixOperationsDevice<Policy, Allocator>::axpy( AMP::Scalar alpha_in,
                                                         const MatrixData &X,
                                                         MatrixData &Y )
{
    auto csrDataX = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( X ) );
    auto csrDataY = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( Y ) );

    AMP_DEBUG_ASSERT( csrDataX );
    AMP_DEBUG_ASSERT( csrDataY );

    AMP_DEBUG_INSIST( csrDataX->d_memory_location != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsDevice is not implemented for device memory" );
    AMP_DEBUG_INSIST( csrDataY->d_memory_location != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsDevice is not implemented for device memory" );
    AMP_DEBUG_INSIST( csrDataX->d_memory_location == csrDataY->d_memory_location,
                      "CSRMatrixOperationsDevice::axpy X and Y must be in same memory space" );

    auto diagMatrixX = csrDataX->getDiagMatrix();
    auto offdMatrixX = csrDataX->getOffdMatrix();

    auto diagMatrixY = csrDataY->getDiagMatrix();
    auto offdMatrixY = csrDataY->getOffdMatrix();

    AMP_DEBUG_ASSERT( diagMatrixX && offdMatrixX );
    AMP_DEBUG_ASSERT( diagMatrixY && offdMatrixY );

    auto alpha = static_cast<scalar_t>( alpha_in );
    CSRLocalMatrixOperationsDevice<Policy, Allocator>::axpy( alpha, diagMatrixX, diagMatrixY );
    if ( csrDataX->hasOffDiag() ) {
        CSRLocalMatrixOperationsDevice<Policy, Allocator>::axpy( alpha, offdMatrixX, offdMatrixY );
    }
}

template<typename Policy, class Allocator>
void CSRMatrixOperationsDevice<Policy, Allocator>::setScalar( AMP::Scalar alpha_in, MatrixData &A )
{
    auto csrData = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );

    AMP_DEBUG_ASSERT( csrData );

    auto diagMatrix = csrData->getDiagMatrix();
    auto offdMatrix = csrData->getOffdMatrix();

    AMP_DEBUG_ASSERT( diagMatrix && offdMatrix );

    AMP_DEBUG_INSIST( csrData->d_memory_location != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsDevice is not implemented for device memory" );

    auto alpha = static_cast<scalar_t>( alpha_in );

    CSRLocalMatrixOperationsDevice<Policy, Allocator>::setScalar( alpha, diagMatrix );
    if ( csrData->hasOffDiag() ) {
        CSRLocalMatrixOperationsDevice<Policy, Allocator>::setScalar( alpha, offdMatrix );
    }
}

template<typename Policy, class Allocator>
void CSRMatrixOperationsDevice<Policy, Allocator>::zero( MatrixData &A )
{
    setScalar( static_cast<scalar_t>( 0.0 ), A );
}

template<typename Policy, class Allocator>
void CSRMatrixOperationsDevice<Policy, Allocator>::setDiagonal( std::shared_ptr<const Vector> in,
                                                                MatrixData &A )
{
    // constrain to one data block for now
    AMP_DEBUG_ASSERT( in && in->numberOfDataBlocks() == 1 && in->isType<scalar_t>( 0 ) );

    const scalar_t *vvals_p = in->getRawDataBlock<scalar_t>();

    auto csrData = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );

    AMP_DEBUG_ASSERT( csrData );

    auto diagMatrix = csrData->getDiagMatrix();

    AMP_DEBUG_ASSERT( diagMatrix );

    AMP_DEBUG_INSIST( csrData->d_memory_location != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsDefault is not implemented for device memory" );

    CSRLocalMatrixOperationsDevice<Policy, Allocator>::setDiagonal( vvals_p, diagMatrix );
}

template<typename Policy, class Allocator>
void CSRMatrixOperationsDevice<Policy, Allocator>::extractDiagonal( MatrixData const &A,
                                                                    std::shared_ptr<Vector> buf )

{
    auto csrData = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );

    AMP_DEBUG_ASSERT( csrData );

    auto diagMatrix = csrData->getDiagMatrix();

    AMP_DEBUG_ASSERT( diagMatrix );

    AMP_DEBUG_INSIST( csrData->d_memory_location != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsDefault is not implemented for device memory" );

    scalar_t *buf_p = buf->getRawDataBlock<scalar_t>();
    CSRLocalMatrixOperationsDevice<Policy, Allocator>::extractDiagonal( diagMatrix, buf_p );
}

template<typename Policy, class Allocator>
void CSRMatrixOperationsDevice<Policy, Allocator>::setIdentity( MatrixData &A )
{
    zero( A );

    auto csrData = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );

    AMP_DEBUG_ASSERT( csrData );

    auto diagMatrix = csrData->getDiagMatrix();

    AMP_DEBUG_ASSERT( diagMatrix );

    AMP_DEBUG_INSIST( csrData->d_memory_location != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsDefault is not implemented for device memory" );

    CSRLocalMatrixOperationsDevice<Policy, Allocator>::setIdentity( diagMatrix );
}

template<typename Policy, class Allocator>
AMP::Scalar CSRMatrixOperationsDevice<Policy, Allocator>::LinfNorm( MatrixData const &A ) const

{
    auto csrData = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );

    AMP_DEBUG_ASSERT( csrData );

    auto diagMatrix = csrData->getDiagMatrix();
    auto offdMatrix = csrData->getOffdMatrix();

    AMP_DEBUG_ASSERT( diagMatrix && offdMatrix );

    AMP_DEBUG_INSIST( csrData->d_memory_location != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsDevice is not implemented for device memory" );

    const auto nRows = csrData->numLocalRows();
    thrust::device_vector<scalar_t> rowSums( nRows, 0.0 );

    CSRLocalMatrixOperationsDevice<Policy, Allocator>::LinfNorm( diagMatrix, rowSums.data().get() );
    if ( csrData->hasOffDiag() ) {
        CSRLocalMatrixOperationsDevice<Policy, Allocator>::LinfNorm( offdMatrix,
                                                                     rowSums.data().get() );
    }

    // Reduce row sums to get global Linf norm
    auto max_norm = *thrust::max_element( thrust::device, rowSums.begin(), rowSums.end() );
    AMP_MPI comm  = csrData->getComm();
    return comm.maxReduce<scalar_t>( max_norm );
}

template<typename Policy, class Allocator>
void CSRMatrixOperationsDevice<Policy, Allocator>::copy( const MatrixData &X, MatrixData &Y )
{
    auto csrDataX = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( X ) );
    auto csrDataY = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( Y ) );

    AMP_DEBUG_ASSERT( csrDataX );
    AMP_DEBUG_ASSERT( csrDataY );

    AMP_DEBUG_INSIST( csrDataX->d_memory_location != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsDevice is not implemented for device memory" );
    AMP_DEBUG_INSIST( csrDataY->d_memory_location != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsDevice is not implemented for device memory" );
    AMP_DEBUG_INSIST( csrDataX->d_memory_location == csrDataY->d_memory_location,
                      "CSRMatrixOperationsDevice::axpy X and Y must be in same memory space" );

    auto diagMatrixX = csrDataX->getDiagMatrix();
    auto offdMatrixX = csrDataX->getOffdMatrix();

    auto diagMatrixY = csrDataY->getDiagMatrix();
    auto offdMatrixY = csrDataY->getOffdMatrix();

    AMP_DEBUG_ASSERT( diagMatrixX && offdMatrixX );
    AMP_DEBUG_ASSERT( diagMatrixY && offdMatrixY );

    CSRLocalMatrixOperationsDevice<Policy, Allocator>::copy( diagMatrixX, diagMatrixY );
    if ( csrDataX->hasOffDiag() ) {
        CSRLocalMatrixOperationsDevice<Policy, Allocator>::copy( offdMatrixX, offdMatrixY );
    }
}

template<typename Policy, class Allocator>
void CSRMatrixOperationsDevice<Policy, Allocator>::copyCast( const MatrixData &X, MatrixData &Y )
{
    auto csrDataY = getCSRMatrixData<Policy, Allocator>( Y );
    AMP_DEBUG_ASSERT( csrDataY );
    if ( X.getCoeffType() == getTypeID<double>() ) {
        using PolicyIn =
            AMP::LinearAlgebra::CSRPolicy<typename Policy::gidx_t, typename Policy::lidx_t, double>;
        auto csrDataX = getCSRMatrixData<PolicyIn, Allocator>( const_cast<MatrixData &>( X ) );
        AMP_DEBUG_ASSERT( csrDataX );

        copyCast<PolicyIn>( csrDataX, csrDataY );
    } else if ( X.getCoeffType() == getTypeID<float>() ) {
        using PolicyIn =
            AMP::LinearAlgebra::CSRPolicy<typename Policy::gidx_t, typename Policy::lidx_t, float>;
        auto csrDataX = getCSRMatrixData<PolicyIn, Allocator>( const_cast<MatrixData &>( X ) );
        AMP_DEBUG_ASSERT( csrDataX );

        copyCast<PolicyIn>( csrDataX, csrDataY );
    } else {
        AMP_ERROR( "Can't copyCast from the given matrix, policy not supported" );
    }
}

template<typename Policy, class Allocator>
template<typename PolicyIn>
void CSRMatrixOperationsDevice<Policy, Allocator>::copyCast( CSRMatrixData<PolicyIn, Allocator> *X,
                                                             matrixdata_t *Y )
{

    AMP_DEBUG_INSIST( X->d_memory_location != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsDevice is not implemented for device memory" );
    AMP_DEBUG_INSIST( Y->d_memory_location != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsDevice is not implemented for device memory" );
    AMP_DEBUG_INSIST( X->d_memory_location == Y->d_memory_location,
                      "CSRMatrixOperationsDevice::copyCast X and Y must be in same memory space" );

    auto diagMatrixX = X->getDiagMatrix();
    auto offdMatrixX = X->getOffdMatrix();

    auto diagMatrixY = Y->getDiagMatrix();
    auto offdMatrixY = Y->getOffdMatrix();

    AMP_DEBUG_ASSERT( diagMatrixX && offdMatrixX );
    AMP_DEBUG_ASSERT( diagMatrixY && offdMatrixY );

    localops_t::template copyCast<PolicyIn>( diagMatrixX, diagMatrixY );
    if ( X->hasOffDiag() ) {
        localops_t::template copyCast<PolicyIn>( offdMatrixX, offdMatrixY );
    }
}

} // namespace AMP::LinearAlgebra

#endif
