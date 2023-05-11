#include "AMP/matrices/DenseSerialMatrix.h"
#include "AMP/matrices/data/DenseSerialMatrixData.h"
#include "AMP/matrices/operations/DenseSerialMatrixOperations.h"
#include "AMP/vectors/VectorBuilder.h"
#include <cstdio>
#include <cstring>

#include <numeric>

namespace AMP::LinearAlgebra {


/********************************************************
 * Constructor/Destructor                                *
 ********************************************************/
DenseSerialMatrix::DenseSerialMatrix( std::shared_ptr<MatrixParameters> params ) : Matrix( params )
{
    d_matrixOps = std::make_shared<DenseSerialMatrixOperations>();
}

DenseSerialMatrix::DenseSerialMatrix( std::shared_ptr<MatrixData> data ) : Matrix( data )
{
    d_matrixOps = std::make_shared<DenseSerialMatrixOperations>();
}

DenseSerialMatrix::~DenseSerialMatrix() {}

/********************************************************
 * Copy/transpose the matrix                             *
 ********************************************************/
std::shared_ptr<Matrix> DenseSerialMatrix::clone() const
{
    return std::make_shared<DenseSerialMatrix>( d_matrixData->cloneMatrixData() );
}

std::shared_ptr<Matrix> DenseSerialMatrix::transpose() const
{
    auto data = d_matrixData->transpose();
    return std::make_shared<DenseSerialMatrix>( data );
}

/********************************************************
 * Multiply two matricies                                *
 * result = this * other_op                              *
 * C(N,M) = A(N,K)*B(K,M)
 ********************************************************/
void DenseSerialMatrix::multiply( std::shared_ptr<Matrix> other_op,
                                  std::shared_ptr<Matrix> &result )
{
    // Create the matrix
    auto params = std::make_shared<AMP::LinearAlgebra::MatrixParameters>(
        getLeftDOFManager(), other_op->getRightDOFManager(), getComm() );

    auto data = std::dynamic_pointer_cast<DenseSerialMatrixData>( d_matrixData );
    AMP_ASSERT( data );
    params->d_VariableLeft  = data->getLeftVariable();
    params->d_VariableRight = data->getRightVariable();
    // Create the matrix
    auto newData   = std::make_shared<AMP::LinearAlgebra::DenseSerialMatrixData>( params );
    auto newMatrix = std::make_shared<AMP::LinearAlgebra::DenseSerialMatrix>( newData );
    AMP_ASSERT( newMatrix );
    result = newMatrix;

    d_matrixOps->matMultiply( *getMatrixData(), *other_op->getMatrixData(), *newData );
}

/********************************************************
 * Get/Set the diagonal                                  *
 ********************************************************/
Vector::shared_ptr DenseSerialMatrix::extractDiagonal( Vector::shared_ptr buf ) const
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
Vector::shared_ptr DenseSerialMatrix::getRightVector() const
{
    auto var = std::dynamic_pointer_cast<DenseSerialMatrixData>( d_matrixData )->getRightVariable();
    return createVector( getRightDOFManager(), var );
}
Vector::shared_ptr DenseSerialMatrix::getLeftVector() const
{
    auto var = std::dynamic_pointer_cast<DenseSerialMatrixData>( d_matrixData )->getLeftVariable();
    return createVector( getLeftDOFManager(), var );
}


} // namespace AMP::LinearAlgebra
