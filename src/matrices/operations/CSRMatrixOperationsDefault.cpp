#include "AMP/matrices/operations/CSRMatrixOperationsDefault.h"
#include "AMP/matrices/data/CSRMatrixData.h"

namespace AMP::LinearAlgebra {
#if 0
static CSRMatrixData const *getCSRMatrixData( MatrixData const &A )
{
    auto ptr = dynamic_cast<CSRMatrixData const *>( &A );
    AMP_INSIST( ptr, "dynamic cast from const MatrixData to const CSRMatrixData failed" );
    return ptr;
}

static CSRMatrixData *getCSRMatrixData( MatrixData &A )
{
    auto ptr = dynamic_cast<CSRMatrixData *>( &A );
    AMP_INSIST( ptr, "dynamic cast from const MatrixData to const CSRMatrixData failed" );
    return ptr;
}
#endif

void CSRMatrixOperationsDefault::mult( std::shared_ptr<const Vector> in,
                                       MatrixData const &A,
                                       std::shared_ptr<Vector> out )
{
    AMP_ERROR( "Not implemented" );
}

void CSRMatrixOperationsDefault::multTranspose( std::shared_ptr<const Vector> in,
                                                MatrixData const &A,
                                                std::shared_ptr<Vector> out )
{
    AMP_ERROR( "Not implemented" );
}

void CSRMatrixOperationsDefault::scale( AMP::Scalar alpha_in, MatrixData &A )
{
    AMP_ERROR( "Not implemented" );
}

void CSRMatrixOperationsDefault::matMultiply( MatrixData const &Am,
                                              MatrixData const &Bm,
                                              MatrixData &Cm )
{
    AMP_ERROR( "Not implemented" );
}

void CSRMatrixOperationsDefault::axpy( AMP::Scalar alpha_in, const MatrixData &X, MatrixData &Y )
{
    AMP_ERROR( "Not implemented" );
}

void CSRMatrixOperationsDefault::setScalar( AMP::Scalar alpha_in, MatrixData &A )
{
    AMP_ERROR( "Not implemented" );
}

void CSRMatrixOperationsDefault::zero( MatrixData &A ) { AMP_ERROR( "Not implemented" ); }

void CSRMatrixOperationsDefault::setDiagonal( std::shared_ptr<const Vector> in, MatrixData &A )
{
    AMP_ERROR( "Not implemented" );
}

void CSRMatrixOperationsDefault::setIdentity( MatrixData &A ) { AMP_ERROR( "Not implemented" ); }

AMP::Scalar CSRMatrixOperationsDefault::L1Norm( MatrixData const &A ) const
{
    AMP_ERROR( "Not implemented" );
}

} // namespace AMP::LinearAlgebra
