#ifndef included_AMP_CSRMatrix_hpp
#define included_AMP_CSRMatrix_hpp

#include "AMP/matrices/CSRMatrix.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/CSRMatrixOperationsDefault.h"
#include "AMP/vectors/VectorBuilder.h"
#include <cstdio>
#include <cstring>

#include <numeric>

namespace AMP::LinearAlgebra {


/********************************************************
 * Constructor/Destructor                                *
 ********************************************************/
template<typename Policy, typename Allocator>
CSRMatrix<Policy,Allocator>::CSRMatrix( std::shared_ptr<MatrixParametersBase> params,
					const std::function<std::vector<size_t>( size_t )> getRow )
  : Matrix( params )
{
    d_matrixOps  = std::make_shared<CSRMatrixOperationsDefault<Policy,Allocator>>();
    d_matrixData = std::make_shared<CSRMatrixData<Policy,Allocator>>( params, getRow );
}

template<typename Policy, typename Allocator>
CSRMatrix<Policy,Allocator>::CSRMatrix( std::shared_ptr<MatrixData> data )
  : Matrix( data )
{
    d_matrixOps = std::make_shared<CSRMatrixOperationsDefault<Policy,Allocator>>();
}

template<typename Policy, typename Allocator>
CSRMatrix<Policy,Allocator>::~CSRMatrix()
{
}

/********************************************************
 * Copy/transpose the matrix                             *
 ********************************************************/
template<typename Policy, typename Allocator>
std::shared_ptr<Matrix> CSRMatrix<Policy,Allocator>::clone() const
{
    return std::make_shared<CSRMatrix<Policy,Allocator>>( d_matrixData->cloneMatrixData() );
}

template<typename Policy, typename Allocator>
std::shared_ptr<Matrix> CSRMatrix<Policy,Allocator>::transpose() const
{
    auto data = d_matrixData->transpose();
    return std::make_shared<CSRMatrix<Policy,Allocator>>( data );
}

/********************************************************
 * Multiply two matricies                                *
 * result = this * other_op                              *
 * C(N,M) = A(N,K)*B(K,M)
 ********************************************************/
template<typename Policy, typename Allocator>
void CSRMatrix<Policy,Allocator>::multiply( std::shared_ptr<Matrix>,
					    std::shared_ptr<Matrix> & )
{
  AMP_ERROR( "Not implemented" );
}

/********************************************************
 * Get/Set the diagonal                                  *
 ********************************************************/
template<typename Policy, typename Allocator>
Vector::shared_ptr CSRMatrix<Policy,Allocator>::extractDiagonal( Vector::shared_ptr buf ) const
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
Vector::shared_ptr CSRMatrix<Policy,Allocator>::getRightVector() const
{
    auto var = std::dynamic_pointer_cast<CSRMatrixData<Policy,Allocator>>( d_matrixData )->getRightVariable();
    return createVector( getRightDOFManager(), var );
}
template<typename Policy, typename Allocator>
Vector::shared_ptr CSRMatrix<Policy,Allocator>::getLeftVector() const
{
    auto var = std::dynamic_pointer_cast<CSRMatrixData<Policy,Allocator>>( d_matrixData )->getLeftVariable();
    return createVector( getLeftDOFManager(), var );
}


} // namespace AMP::LinearAlgebra

#endif
