#include "AMP/matrices/DenseSerialMatrix.h"
#include "AMP/matrices/data/DenseSerialMatrixData.h"
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
}

DenseSerialMatrix::DenseSerialMatrix( std::shared_ptr<MatrixData> data ) : Matrix( data ) {}

DenseSerialMatrix::~DenseSerialMatrix() {}

/********************************************************
 * Copy/transpose the matrix                             *
 ********************************************************/
std::shared_ptr<Matrix> DenseSerialMatrix::cloneMatrix() const
{
    return std::make_shared<DenseSerialMatrix>( d_matrixData->cloneMatrixData() );
}

std::shared_ptr<Matrix> DenseSerialMatrix::transpose() const
{
    auto M2 = cloneMatrix();

    auto m2Data       = std::dynamic_pointer_cast<DenseSerialMatrixData>( M2->getMatrixData() );
    const auto m1Data = std::dynamic_pointer_cast<const DenseSerialMatrixData>( getMatrixData() );
    AMP_ASSERT( m1Data && m2Data );
    const auto nrows      = m1Data->d_rows;
    const auto ncols      = m1Data->d_cols;
    auto *m2RawData       = m2Data->d_M;
    const auto *m1RawData = m1Data->d_M;

    for ( size_t i = 0; i < nrows; ++i ) {
        for ( size_t j = 0; j < ncols; ++j ) {
            m2RawData[j + i * ncols] = m1RawData[i + j * nrows];
        }
    }
    return M2;
}


/********************************************************
 * Matrix-vector multiplication                          *
 ********************************************************/
void DenseSerialMatrix::mult( AMP::LinearAlgebra::Vector::const_shared_ptr in,
                              AMP::LinearAlgebra::Vector::shared_ptr out )
{
    auto m1Data = std::dynamic_pointer_cast<DenseSerialMatrixData>( getMatrixData() );
    AMP_ASSERT( m1Data );
    const auto nrows      = m1Data->d_rows;
    const auto ncols      = m1Data->d_cols;
    const auto *m1RawData = m1Data->d_M;

    AMP_ASSERT( in->getGlobalSize() == ncols );
    AMP_ASSERT( out->getGlobalSize() == nrows );
    auto k = new size_t[std::max( ncols, nrows )];
    for ( size_t i = 0; i < std::max( ncols, nrows ); i++ )
        k[i] = i;
    // Get x
    auto x = new double[ncols];
    in->getValuesByGlobalID( ncols, k, x );
    // Initialize y
    auto y = new double[nrows];
    memset( y, 0, nrows * sizeof( double ) );
    // Perform y = M*x
    for ( size_t j = 0; j < ncols; j++ ) {
        for ( size_t i = 0; i < nrows; i++ )
            y[i] += m1RawData[i + j * nrows] * x[j];
    }
    // Save y
    out->setValuesByGlobalID( nrows, k, y );
    delete[] x;
    delete[] y;
    delete[] k;
}
void DenseSerialMatrix::multTranspose( AMP::LinearAlgebra::Vector::const_shared_ptr in,
                                       AMP::LinearAlgebra::Vector::shared_ptr out )
{
    auto m1Data = std::dynamic_pointer_cast<DenseSerialMatrixData>( getMatrixData() );
    AMP_ASSERT( m1Data );
    const auto nrows      = m1Data->d_rows;
    const auto ncols      = m1Data->d_cols;
    const auto *m1RawData = m1Data->d_M;

    AMP_ASSERT( in->getGlobalSize() == nrows );
    AMP_ASSERT( out->getGlobalSize() == ncols );
    auto k = new size_t[std::max( ncols, nrows )];
    for ( size_t i = 0; i < std::max( ncols, nrows ); i++ )
        k[i] = i;
    // Get x
    auto x = new double[nrows];
    in->getValuesByGlobalID( nrows, k, x );
    // Initialize y
    auto y = new double[ncols];
    memset( y, 0, ncols * sizeof( double ) );
    // Perform y = M*x
    for ( size_t j = 0; j < ncols; j++ ) {
        for ( size_t i = 0; i < nrows; i++ )
            y[j] += m1RawData[i + j * nrows] * x[i];
    }
    // Save y
    out->setValuesByGlobalID( ncols, k, y );
    delete[] x;
    delete[] y;
    delete[] k;
}


/********************************************************
 * Scale/axpy/setScalar                                  *
 ********************************************************/
void DenseSerialMatrix::scale( double alpha )
{
    auto m1Data = std::dynamic_pointer_cast<DenseSerialMatrixData>( getMatrixData() );
    AMP_ASSERT( m1Data );
    const auto nrows = m1Data->d_rows;
    const auto ncols = m1Data->d_cols;
    auto *m1RawData  = m1Data->d_M;

    for ( size_t i = 0; i < nrows * ncols; i++ )
        m1RawData[i] *= alpha;
}

void DenseSerialMatrix::axpy( double alpha, const Matrix &X )
{
    auto m1Data = std::dynamic_pointer_cast<DenseSerialMatrixData>( getMatrixData() );
    AMP_ASSERT( m1Data );
    auto *m1RawData  = m1Data->d_M;
    const auto nrows = m1Data->d_rows;
    const auto ncols = m1Data->d_cols;

    AMP_ASSERT( X.numGlobalRows() == this->numGlobalRows() );
    AMP_ASSERT( X.numGlobalColumns() == this->numGlobalColumns() );
    if ( dynamic_cast<const DenseSerialMatrix *>( &X ) == nullptr ) {
        // X is an unknown matrix type
        std::vector<size_t> cols;
        std::vector<double> values;
        for ( size_t i = 0; i < nrows; i++ ) {
            X.getRowByGlobalID( static_cast<int>( i ), cols, values );
            for ( size_t j = 0; j < cols.size(); j++ )
                m1RawData[i + cols[j] * nrows] += alpha * values[j];
        }
    } else {
        // We are dealing with two DenseSerialMatrix classes
        const auto m2Data =
            std::dynamic_pointer_cast<const DenseSerialMatrixData>( X.getMatrixData() );
        AMP_ASSERT( m2Data );
        auto *m2RawData = m2Data->d_M;
        for ( size_t i = 0; i < nrows * ncols; i++ ) {
            m1RawData[i] += alpha * m2RawData[i];
        }
    }
}
void DenseSerialMatrix::setScalar( double alpha )
{
    auto m1Data = std::dynamic_pointer_cast<DenseSerialMatrixData>( getMatrixData() );
    AMP_ASSERT( m1Data );
    const auto nrows = m1Data->d_rows;
    const auto ncols = m1Data->d_cols;
    auto *m1RawData  = m1Data->d_M;

    for ( size_t i = 0; i < nrows * ncols; i++ )
        m1RawData[i] = alpha;
}
void DenseSerialMatrix::zero()
{
    auto m1Data = std::dynamic_pointer_cast<DenseSerialMatrixData>( getMatrixData() );
    AMP_ASSERT( m1Data );
    const auto nrows = m1Data->d_rows;
    const auto ncols = m1Data->d_cols;
    auto *m1RawData  = m1Data->d_M;

    memset( m1RawData, 0, nrows * ncols * sizeof( double ) );
}

/********************************************************
 * Get/Set the diagonal                                  *
 ********************************************************/
Vector::shared_ptr DenseSerialMatrix::extractDiagonal( Vector::shared_ptr buf ) const
{
    const auto m1Data = std::dynamic_pointer_cast<const DenseSerialMatrixData>( getMatrixData() );
    AMP_ASSERT( m1Data );
    const auto nrows      = m1Data->d_rows;
    const auto ncols      = m1Data->d_cols;
    const auto *m1RawData = m1Data->d_M;

    AMP_ASSERT( ncols == nrows );
    Vector::shared_ptr out = buf;
    if ( buf == nullptr )
        out = this->getRightVector();
    AMP_ASSERT( out->getGlobalSize() == ncols );
    auto y = new double[ncols];
    for ( size_t i = 0; i < ncols; i++ )
        y[i] = m1RawData[i + i * nrows];
    auto k = new size_t[ncols];
    for ( size_t i = 0; i < ncols; i++ )
        k[i] = i;
    out->setValuesByGlobalID( ncols, k, y );
    delete[] y;
    delete[] k;
    return out;
}
void DenseSerialMatrix::setDiagonal( Vector::const_shared_ptr in )
{
    auto m1Data = std::dynamic_pointer_cast<DenseSerialMatrixData>( getMatrixData() );
    AMP_ASSERT( m1Data );
    const auto nrows = m1Data->d_rows;
    const auto ncols = m1Data->d_cols;
    auto *m1RawData  = m1Data->d_M;

    AMP_ASSERT( ncols == nrows );
    AMP_ASSERT( in->getGlobalSize() == nrows );
    auto k = new size_t[nrows];
    for ( size_t i = 0; i < nrows; i++ )
        k[i] = i;
    auto x = new double[nrows];
    in->getValuesByGlobalID( nrows, k, x );
    for ( size_t i = 0; i < nrows; i++ )
        m1RawData[i + i * nrows] = x[i];
    delete[] x;
}
void DenseSerialMatrix::setIdentity()
{
    auto m1Data = std::dynamic_pointer_cast<DenseSerialMatrixData>( getMatrixData() );
    AMP_ASSERT( m1Data );
    const auto nrows = m1Data->d_rows;
    const auto ncols = m1Data->d_cols;
    auto *m1RawData  = m1Data->d_M;

    AMP_ASSERT( ncols == nrows );
    memset( m1RawData, 0, nrows * ncols * sizeof( double ) );
    for ( size_t i = 0; i < nrows; i++ )
        m1RawData[i + i * nrows] = 1.0;
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


/********************************************************
 * Compute the maximum column sum                        *
 ********************************************************/
double DenseSerialMatrix::L1Norm() const
{
    auto m1Data = std::dynamic_pointer_cast<const DenseSerialMatrixData>( getMatrixData() );
    AMP_ASSERT( m1Data );
    const auto nrows = m1Data->d_rows;
    const auto ncols = m1Data->d_cols;
    auto *m1RawData  = m1Data->d_M;

    double norm = 0.0;
    for ( size_t j = 0; j < ncols; j++ ) {
        double sum = 0.0;
        for ( size_t i = 0; i < nrows; i++ )
            sum += fabs( m1RawData[i + j * nrows] );
        norm = std::max( norm, sum );
    }
    return norm;
}

/********************************************************
 * Multiply two matricies                                *
 * result = this * other_op                              *
 * C(N,M) = A(N,K)*B(K,M)
 ********************************************************/
void DenseSerialMatrix::multiply( std::shared_ptr<Matrix> other_op,
                                  std::shared_ptr<Matrix> &result )
{
    if ( this->numGlobalColumns() != other_op->numGlobalRows() )
        AMP_ERROR( "Inner matrix dimensions must agree" );
    size_t N = this->numGlobalRows();
    size_t K = this->numGlobalColumns();
    size_t M = other_op->numGlobalColumns();
    // Create the matrix
    auto params = std::make_shared<AMP::LinearAlgebra::MatrixParameters>(
        other_op->getRightDOFManager(), getRightDOFManager(), getComm() );

    auto data = std::dynamic_pointer_cast<DenseSerialMatrixData>( d_matrixData );
    AMP_ASSERT( data );
    params->d_VariableLeft  = data->getLeftVariable();
    params->d_VariableRight = data->getRightVariable();
    // Create the matrix
    auto newData   = std::make_shared<AMP::LinearAlgebra::DenseSerialMatrixData>( params );
    auto newMatrix = std::make_shared<AMP::LinearAlgebra::DenseSerialMatrix>( newData );
    AMP_ASSERT( newMatrix );
    result    = newMatrix;
    double *C = std::dynamic_pointer_cast<DenseSerialMatrixData>( newMatrix->getMatrixData() )->d_M;
    memset( C, 0, N * M * sizeof( double ) );
    // Perform the muliplication
    if ( std::dynamic_pointer_cast<DenseSerialMatrix>( other_op ) == nullptr ) {
        // X is an unknown matrix type
        AMP_ERROR( "Not programmed yet" );
    } else {
        // We are dealing with all DenseSerialMatrix classes
        const double *A = data->d_M;
        const double *B =
            std::dynamic_pointer_cast<DenseSerialMatrixData>( other_op->getMatrixData() )->d_M;
        for ( size_t m = 0; m < M; m++ ) {
            for ( size_t k = 0; k < K; k++ ) {
                double b = B[k + m * K];
                for ( size_t n = 0; n < N; n++ ) {
                    C[n + m * N] += A[n + k * N] * b;
                }
            }
        }
    }
}


} // namespace AMP::LinearAlgebra
