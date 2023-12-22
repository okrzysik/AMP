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
CSRMatrix::CSRMatrix( std::shared_ptr<MatrixParameters> params ) : Matrix( params )
{
    d_matrixOps  = std::make_shared<CSRMatrixOperationsDefault>();
    d_matrixData = std::make_shared<CSRMatrixData>( params );
}

CSRMatrix::CSRMatrix( std::shared_ptr<MatrixData> data ) : Matrix( data )
{
    d_matrixOps = std::make_shared<CSRMatrixOperationsDefault>();
}

CSRMatrix::~CSRMatrix() {}

/********************************************************
 * Copy/transpose the matrix                             *
 ********************************************************/
std::shared_ptr<Matrix> CSRMatrix::clone() const
{
    return std::make_shared<CSRMatrix>( d_matrixData->cloneMatrixData() );
}

std::shared_ptr<Matrix> CSRMatrix::transpose() const
{
    auto data = d_matrixData->transpose();
    return std::make_shared<CSRMatrix>( data );
}

/********************************************************
 * Multiply two matricies                                *
 * result = this * other_op                              *
 * C(N,M) = A(N,K)*B(K,M)
 ********************************************************/
void CSRMatrix::multiply( std::shared_ptr<Matrix> other_op, std::shared_ptr<Matrix> &result )
{
    // Create the matrix
    auto params = std::make_shared<AMP::LinearAlgebra::MatrixParameters>(
        getLeftDOFManager(), other_op->getRightDOFManager(), getComm() );

    // Create the matrix
    auto newData   = std::make_shared<AMP::LinearAlgebra::CSRMatrixData>( params );
    auto newMatrix = std::make_shared<AMP::LinearAlgebra::CSRMatrix>( newData );
    AMP_ASSERT( newMatrix );
    result = newMatrix;

    d_matrixOps->matMultiply( *getMatrixData(), *other_op->getMatrixData(), *newData );
}

/********************************************************
 * Get/Set the diagonal                                  *
 ********************************************************/
Vector::shared_ptr CSRMatrix::extractDiagonal( Vector::shared_ptr buf ) const
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
Vector::shared_ptr CSRMatrix::getRightVector() const { AMP_ERROR( "Not implemented" ); }
Vector::shared_ptr CSRMatrix::getLeftVector() const { AMP_ERROR( "Not implemented" ); }


} // namespace AMP::LinearAlgebra
