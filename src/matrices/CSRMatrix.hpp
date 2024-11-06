#ifndef included_AMP_CSRMatrix_hpp
#define included_AMP_CSRMatrix_hpp

#include "AMP/matrices/CSRMatrix.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/default/CSRMatrixOperationsDefault.h"
#if defined( USE_DEVICE ) && ( defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS ) )
    #include "AMP/matrices/operations/kokkos/CSRMatrixOperationsKokkos.h"
#endif
#include "AMP/utils/memory.h"
#include "AMP/vectors/VectorBuilder.h"

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
#if defined( USE_DEVICE ) && ( defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS ) )
    if ( std::is_same<Allocator, AMP::HostAllocator<int>>::value ) {
        AMP::pout << "Using default operations" << std::endl;
        d_matrixOps = std::make_shared<CSRMatrixOperationsDefault<Policy, Allocator>>();
    } else {
        AMP::pout << "Using Kokkos operations" << std::endl;
        d_matrixOps = std::make_shared<CSRMatrixOperationsKokkos<Policy, Allocator>>();
    }
#else
    AMP::pout << "Using default operations" << std::endl;
    d_matrixOps = std::make_shared<CSRMatrixOperationsDefault<Policy, Allocator>>();
#endif
    d_matrixData = std::make_shared<CSRMatrixData<Policy, Allocator>>( params );
}

template<typename Policy, typename Allocator>
CSRMatrix<Policy, Allocator>::CSRMatrix( std::shared_ptr<MatrixData> data ) : Matrix( data )
{
#if defined( USE_DEVICE ) && ( defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS ) )
    if ( std::is_same<Allocator, AMP::HostAllocator<int>>::value ) {
        AMP::pout << "Using default operations" << std::endl;
        d_matrixOps = std::make_shared<CSRMatrixOperationsDefault<Policy, Allocator>>();
    } else {
        AMP::pout << "Using Kokkos operations" << std::endl;
        d_matrixOps = std::make_shared<CSRMatrixOperationsKokkos<Policy, Allocator>>();
    }
#else
    AMP::pout << "Using default operations" << std::endl;
    d_matrixOps = std::make_shared<CSRMatrixOperationsDefault<Policy, Allocator>>();
#endif
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
    auto data = d_matrixData->transpose();
    return std::make_shared<CSRMatrix<Policy, Allocator>>( data );
}

/********************************************************
 * Multiply two matricies                                *
 * result = this * other_op                              *
 * C(N,M) = A(N,K)*B(K,M)
 ********************************************************/
template<typename Policy, typename Allocator>
void CSRMatrix<Policy, Allocator>::multiply( std::shared_ptr<Matrix> other_op,
                                             std::shared_ptr<Matrix> &result )
{
    // Create the matrix
    auto params = std::make_shared<AMP::LinearAlgebra::MatrixParameters>(
        getLeftDOFManager(), other_op->getRightDOFManager(), getComm() );

    // Create the matrix
    auto newData = std::make_shared<AMP::LinearAlgebra::CSRMatrixData<Policy, Allocator>>( params );
    auto newMatrix = std::make_shared<AMP::LinearAlgebra::CSRMatrix<Policy, Allocator>>( newData );
    AMP_ASSERT( newMatrix );
    result = newMatrix;

    d_matrixOps->matMultiply( *getMatrixData(), *other_op->getMatrixData(), *newData );
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
