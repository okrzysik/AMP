#ifndef included_AMP_CSRMatrix_hpp
#define included_AMP_CSRMatrix_hpp

#include "AMP/vectors/VectorBuilder.h"
#include "AMP/matrices/CSRMatrix.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/CSRMatrixOperationsDefault.h"

#if defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS )
#include "Kokkos_Core.hpp"
#include "AMP/matrices/operations/CSRMatrixOperationsKokkos.h"
#endif

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
#if defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS )
    d_matrixOps  = std::make_shared<CSRMatrixOperationsKokkos<Policy, Allocator, Kokkos::DefaultExecutionSpace>>();
#else
    d_matrixOps  = std::make_shared<CSRMatrixOperationsDefault<Policy, Allocator>>();
#endif
    d_matrixData = std::make_shared<CSRMatrixData<Policy, Allocator>>( params );
}

template<typename Policy, typename Allocator>
CSRMatrix<Policy, Allocator>::CSRMatrix( std::shared_ptr<MatrixData> data ) : Matrix( data )
{
#if defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS )
    d_matrixOps  = std::make_shared<CSRMatrixOperationsKokkos<Policy, Allocator, Kokkos::DefaultExecutionSpace>>();
#else
  d_matrixOps  = std::make_shared<CSRMatrixOperationsDefault<Policy, Allocator>>();
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
void CSRMatrix<Policy, Allocator>::multiply( std::shared_ptr<Matrix>, std::shared_ptr<Matrix> & )
{
    AMP_ERROR( "Not implemented" );
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

    d_matrixData->extractDiagonal( out );

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
    return createVector( getRightDOFManager(), var );
}
template<typename Policy, typename Allocator>
Vector::shared_ptr CSRMatrix<Policy, Allocator>::getLeftVector() const
{
    auto var = std::dynamic_pointer_cast<CSRMatrixData<Policy, Allocator>>( d_matrixData )
                   ->getLeftVariable();
    return createVector( getLeftDOFManager(), var );
}


} // namespace AMP::LinearAlgebra

#endif
