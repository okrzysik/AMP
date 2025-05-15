#ifndef included_AMP_CSRMatrix_hpp
#define included_AMP_CSRMatrix_hpp

#include "AMP/matrices/CSRMatrix.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/default/CSRMatrixOperationsDefault.h"

#ifdef USE_DEVICE
    #include "AMP/matrices/operations/device/CSRMatrixOperationsDevice.h"
#endif

#if ( defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS ) )
    #include "AMP/matrices/operations/kokkos/CSRMatrixOperationsKokkos.h"
#endif

#include "AMP/utils/memory.h"
#include "AMP/vectors/VectorBuilder.h"

#include "ProfilerApp.h"

#include <cstdio>
#include <cstring>
#include <numeric>

namespace AMP::LinearAlgebra {

/********************************************************
 * Constructor/Destructor                                *
 ********************************************************/
template<typename Policy, typename Allocator>
CSRMatrix<Policy, Allocator>::CSRMatrix( std::shared_ptr<MatrixParametersBase> params )
    : Matrix( params )
{
    bool set_ops = false;

    if ( params->d_backend == AMP::Utilities::Backend::hip_cuda ) {
#ifdef USE_DEVICE
        d_matrixOps = std::make_shared<CSRMatrixOperationsDevice<Policy, Allocator>>();
        set_ops     = true;
#else
        AMP_ERROR( "HIP or CUDA need to be loaded to be able to use hip_cuda backend." );
#endif
    }

    if ( params->d_backend == AMP::Utilities::Backend::kokkos ) {
#if ( defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS ) )
        d_matrixOps = std::make_shared<CSRMatrixOperationsKokkos<Policy, Allocator>>();
        set_ops     = true;
#else
        AMP_ERROR( "KOKKOS need to be loaded to be able to use kokkos backend." );
#endif
    }

    // nothing above matched, fall back on default operations
    if ( !set_ops ) {
        d_matrixOps = std::make_shared<CSRMatrixOperationsDefault<Policy, Allocator>>();
    }

    d_matrixData = std::make_shared<CSRMatrixData<Policy, Allocator>>( params );
}

template<typename Policy, typename Allocator>
CSRMatrix<Policy, Allocator>::CSRMatrix( std::shared_ptr<MatrixData> data ) : Matrix( data )
{
    auto backend = data->getBackend();
    bool set_ops = false;

    if ( backend == AMP::Utilities::Backend::hip_cuda ) {
#ifdef USE_DEVICE
        d_matrixOps = std::make_shared<CSRMatrixOperationsDevice<Policy, Allocator>>();
        set_ops     = true;
#else
        AMP_ERROR( "HIP or CUDA need to be loaded to be able to use hip_cuda backend." );
#endif
    }

    if ( backend == AMP::Utilities::Backend::kokkos ) {
#if ( defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS ) )
        d_matrixOps = std::make_shared<CSRMatrixOperationsKokkos<Policy, Allocator>>();
        set_ops     = true;
#else
        AMP_ERROR( "KOKKOS need to be loaded to be able to use kokkos backend." );
#endif
    }

    // nothing above matched, fall back on default operations
    if ( !set_ops ) {
        d_matrixOps = std::make_shared<CSRMatrixOperationsDefault<Policy, Allocator>>();
    }
}

template<typename Policy, typename Allocator>
CSRMatrix<Policy, Allocator>::~CSRMatrix()
{
}

/********************************************************
 * Copy/transpose the matrix                             *
 ********************************************************/
template<typename Policy, typename Allocator>
std::shared_ptr<Matrix> CSRMatrix<Policy, Allocator>::clone() const
{
    return std::make_shared<CSRMatrix<Policy, Allocator>>( d_matrixData->cloneMatrixData() );
}

template<typename Policy, typename Allocator>
std::shared_ptr<Matrix> CSRMatrix<Policy, Allocator>::transpose() const
{
    PROFILE( "CSRMatrix<Policy, Allocator>::transpose" );

    auto data = d_matrixData->transpose();
    return std::make_shared<CSRMatrix<Policy, Allocator>>( data );
}

/********************************************************
 * Multiply two matrices                                *
 * result = this * other_op                             *
 * C(N,M) = A(N,K)*B(K,M)
 ********************************************************/
template<typename Policy, typename Allocator>
void CSRMatrix<Policy, Allocator>::multiply( std::shared_ptr<Matrix> other_op,
                                             std::shared_ptr<Matrix> &result )
{
    PROFILE( "CSRMatrix<Policy, Allocator>::multiply" );

    // if the result is empty then create it
    if ( result.get() == nullptr ) {
        // pull out matrix data objects and ensure they are of correct type
        auto thisData = std::dynamic_pointer_cast<CSRMatrixData<Policy, Allocator>>( d_matrixData );
        auto otherData = std::dynamic_pointer_cast<CSRMatrixData<Policy, Allocator>>(
            other_op->getMatrixData() );
        AMP_DEBUG_INSIST( thisData && otherData,
                          "CSRMatrix::multiply received invalid MatrixData types" );

        // Build matrix parameters object for result from this op and the other op
        auto params = std::make_shared<AMP::LinearAlgebra::MatrixParameters>(
            getLeftDOFManager(),
            other_op->getRightDOFManager(),
            getComm(),
            thisData->getLeftVariable(),
            otherData->getRightVariable(),
            std::function<std::vector<size_t>( size_t )>() );

        // Create the matrix
        auto newData =
            std::make_shared<AMP::LinearAlgebra::CSRMatrixData<Policy, Allocator>>( params );
        std::shared_ptr<Matrix> newMatrix =
            std::make_shared<AMP::LinearAlgebra::CSRMatrix<Policy, Allocator>>( newData );
        AMP_ASSERT( newMatrix );
        result.swap( newMatrix );
    } else {
        // do something to check that result is compatible with this and other_op?
    }

    d_matrixOps->matMultiply(
        *getMatrixData(), *other_op->getMatrixData(), *result->getMatrixData() );
}

/********************************************************
 * Get/Set the diagonal                                  *
 ********************************************************/
template<typename Policy, typename Allocator>
Vector::shared_ptr CSRMatrix<Policy, Allocator>::extractDiagonal( Vector::shared_ptr buf ) const
{
    Vector::shared_ptr out = buf;
    if ( !buf )
        out = this->getRightVector();

    d_matrixOps->extractDiagonal( *getMatrixData(), out );

    return out;
}

/********************************************************
 * Get the left/right vectors and DOFManagers            *
 ********************************************************/
template<typename Policy, typename Allocator>
Vector::shared_ptr CSRMatrix<Policy, Allocator>::getRightVector() const
{
    auto var = std::dynamic_pointer_cast<CSRMatrixData<Policy, Allocator>>( d_matrixData )
                   ->getRightVariable();
    const auto memloc = AMP::Utilities::getAllocatorMemoryType<Allocator>();
    return createVector( getRightDOFManager(), var, true, memloc );
}

template<typename Policy, typename Allocator>
Vector::shared_ptr CSRMatrix<Policy, Allocator>::getLeftVector() const
{
    auto var = std::dynamic_pointer_cast<CSRMatrixData<Policy, Allocator>>( d_matrixData )
                   ->getLeftVariable();
    const auto memloc = AMP::Utilities::getAllocatorMemoryType<Allocator>();
    return createVector( getLeftDOFManager(), var, true, memloc );
}


} // namespace AMP::LinearAlgebra

#endif
