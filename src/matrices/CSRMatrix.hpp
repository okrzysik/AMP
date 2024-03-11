#ifndef included_AMP_CSRMatrix_hpp
#define included_AMP_CSRMatrix_hpp

#include "AMP/matrices/CSRMatrix.h"
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
template<typename Policy>
CSRMatrix<Policy>::CSRMatrix( std::shared_ptr<MatrixParametersBase> params ) : Matrix( params )
{
    d_matrixOps  = std::make_shared<CSRMatrixOperationsDefault<Policy>>();
    d_matrixData = std::make_shared<CSRMatrixData<Policy>>( params );
}

template<typename Policy>
CSRMatrix<Policy>::CSRMatrix( std::shared_ptr<MatrixData> data ) : Matrix( data )
{
    d_matrixOps = std::make_shared<CSRMatrixOperationsDefault<Policy>>();
}

template<typename Policy>
CSRMatrix<Policy>::~CSRMatrix()
{
}

/********************************************************
 * Copy/transpose the matrix                             *
 ********************************************************/
template<typename Policy>
std::shared_ptr<Matrix> CSRMatrix<Policy>::clone() const
{
    return std::make_shared<CSRMatrix<Policy>>( d_matrixData->cloneMatrixData() );
}

template<typename Policy>
std::shared_ptr<Matrix> CSRMatrix<Policy>::transpose() const
{
    auto data = d_matrixData->transpose();
    return std::make_shared<CSRMatrix<Policy>>( data );
}

/********************************************************
 * Multiply two matricies                                *
 * result = this * other_op                              *
 * C(N,M) = A(N,K)*B(K,M)
 ********************************************************/
template<typename Policy>
void CSRMatrix<Policy>::multiply( std::shared_ptr<Matrix> other_op,
                                  std::shared_ptr<Matrix> &result )
{
    // Create the matrix
    auto params = std::make_shared<AMP::LinearAlgebra::MatrixParameters>(
        getLeftDOFManager(), other_op->getRightDOFManager(), getComm() );

    // Create the matrix
    auto newData   = std::make_shared<AMP::LinearAlgebra::CSRMatrixData<Policy>>( params );
    auto newMatrix = std::make_shared<AMP::LinearAlgebra::CSRMatrix<Policy>>( newData );
    AMP_ASSERT( newMatrix );
    result = newMatrix;

    d_matrixOps->matMultiply( *getMatrixData(), *other_op->getMatrixData(), *newData );
}

/********************************************************
 * Get/Set the diagonal                                  *
 ********************************************************/
template<typename Policy>
Vector::shared_ptr CSRMatrix<Policy>::extractDiagonal( Vector::shared_ptr buf ) const
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
template<typename Policy>
Vector::shared_ptr CSRMatrix<Policy>::getRightVector() const
{
    auto var = std::dynamic_pointer_cast<CSRMatrixData<Policy>>( d_matrixData )->getRightVariable();
    return createVector( getRightDOFManager(), var );
}
template<typename Policy>
Vector::shared_ptr CSRMatrix<Policy>::getLeftVector() const
{
    auto var = std::dynamic_pointer_cast<CSRMatrixData<Policy>>( d_matrixData )->getLeftVariable();
    return createVector( getLeftDOFManager(), var );
}


} // namespace AMP::LinearAlgebra

#endif
