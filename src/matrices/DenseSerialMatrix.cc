#include "matrices/DenseSerialMatrix.h"
#include "string.h"
#include "vectors/VectorBuilder.h"
#include <stdio.h>


namespace AMP {
namespace LinearAlgebra {


/********************************************************
* Constructor/Destructor                                *
********************************************************/
DenseSerialMatrix::DenseSerialMatrix( MatrixParameters::shared_ptr params ) : Matrix( params )
{
    d_VariableLeft    = params->d_VariableLeft;
    d_VariableRight   = params->d_VariableRight;
    d_DOFManagerLeft  = params->getLeftDOFManager();
    d_DOFManagerRight = params->getRightDOFManager();
    d_rows            = params->getGlobalNumberOfRows();
    d_cols            = params->getGlobalNumberOfColumns();
    d_M               = new double[d_rows * d_cols];
    memset( d_M, 0, d_rows * d_cols * sizeof( double ) );
}
DenseSerialMatrix::~DenseSerialMatrix() { delete[] d_M; }


/********************************************************
* Copy/transpose the matrix                             *
********************************************************/
Matrix::shared_ptr DenseSerialMatrix::transpose() const
{
    // Create the matrix parameters
    AMP::shared_ptr<AMP::LinearAlgebra::MatrixParameters> params(
        new AMP::LinearAlgebra::MatrixParameters( d_DOFManagerRight, d_DOFManagerLeft, d_comm ) );
    params->d_VariableLeft  = d_VariableRight;
    params->d_VariableRight = d_VariableLeft;
    // Create the matrix
    AMP::shared_ptr<AMP::LinearAlgebra::DenseSerialMatrix> newMatrix(
        new AMP::LinearAlgebra::DenseSerialMatrix( params ) );
    double *M2 = newMatrix->d_M;
    for ( size_t i = 0; i < d_rows; i++ ) {
        for ( size_t j = 0; j < d_cols; j++ ) {
            M2[j + i * d_cols] = d_M[i + j * d_rows];
        }
    }
    return newMatrix;
}
Matrix::shared_ptr DenseSerialMatrix::cloneMatrix() const
{
    // Create the matrix parameters
    AMP::shared_ptr<AMP::LinearAlgebra::MatrixParameters> params(
        new AMP::LinearAlgebra::MatrixParameters( d_DOFManagerLeft, d_DOFManagerRight, d_comm ) );
    params->d_VariableLeft  = d_VariableLeft;
    params->d_VariableRight = d_VariableRight;
    // Create the matrix
    AMP::shared_ptr<AMP::LinearAlgebra::DenseSerialMatrix> newMatrix(
        new AMP::LinearAlgebra::DenseSerialMatrix( params ) );
    double *M2 = newMatrix->d_M;
    memcpy( M2, d_M, d_cols * d_rows * sizeof( double ) );
    return newMatrix;
}


/********************************************************
* Matrix-vector multiplication                          *
********************************************************/
void DenseSerialMatrix::mult( AMP::LinearAlgebra::Vector::const_shared_ptr in,
                              AMP::LinearAlgebra::Vector::shared_ptr out )
{
    AMP_ASSERT( in->getGlobalSize() == d_cols );
    AMP_ASSERT( out->getGlobalSize() == d_rows );
    size_t *k = new size_t[std::max( d_cols, d_rows )];
    for ( size_t i = 0; i < std::max( d_cols, d_rows ); i++ )
        k[i]       = i;
    // Get x
    double *x = new double[d_cols];
    in->getValuesByGlobalID( d_cols, k, x );
    // Initialize y
    double *y = new double[d_rows];
    memset( y, 0, d_rows * sizeof( double ) );
    // Perform y = M*x
    for ( size_t j = 0; j < d_cols; j++ ) {
        for ( size_t i = 0; i < d_rows; i++ )
            y[i] += d_M[i + j * d_rows] * x[j];
    }
    // Save y
    out->setValuesByGlobalID( d_rows, k, y );
    delete[] x;
    delete[] y;
    delete[] k;
}
void DenseSerialMatrix::multTranspose( AMP::LinearAlgebra::Vector::const_shared_ptr in,
                                       AMP::LinearAlgebra::Vector::shared_ptr out )
{
    AMP_ASSERT( in->getGlobalSize() == d_rows );
    AMP_ASSERT( out->getGlobalSize() == d_cols );
    size_t *k = new size_t[std::max( d_cols, d_rows )];
    for ( size_t i = 0; i < std::max( d_cols, d_rows ); i++ )
        k[i]       = i;
    // Get x
    double *x = new double[d_rows];
    in->getValuesByGlobalID( d_rows, k, x );
    // Initialize y
    double *y = new double[d_cols];
    memset( y, 0, d_cols * sizeof( double ) );
    // Perform y = M*x
    for ( size_t j = 0; j < d_cols; j++ ) {
        for ( size_t i = 0; i < d_rows; i++ )
            y[j] += d_M[i + j * d_rows] * x[i];
    }
    // Save y
    out->setValuesByGlobalID( d_cols, k, y );
    delete[] x;
    delete[] y;
    delete[] k;
}


/********************************************************
* Scale/axpy/setScalar                                  *
********************************************************/
void DenseSerialMatrix::scale( double alpha )
{
    for ( size_t i = 0; i < d_rows * d_cols; i++ )
        d_M[i] *= alpha;
}
void DenseSerialMatrix::axpy( double alpha, const Matrix &X )
{
    AMP_ASSERT( X.numGlobalRows() == this->numGlobalRows() );
    AMP_ASSERT( X.numGlobalColumns() == this->numGlobalColumns() );
    if ( dynamic_cast<const DenseSerialMatrix *>( &X ) == nullptr ) {
        // X is an unknown matrix type
        std::vector<unsigned int> cols;
        std::vector<double> values;
        for ( size_t i = 0; i < d_rows; i++ ) {
            X.getRowByGlobalID( static_cast<int>( i ), cols, values );
            for ( size_t j = 0; j < cols.size(); j++ )
                d_M[i + cols[j] * d_rows] += alpha * values[j];
        }
    } else {
        // We are dealing with two DenseSerialMatrix classes
        const double *M2 = dynamic_cast<const DenseSerialMatrix *>( &X )->d_M;
        for ( size_t i = 0; i < d_rows * d_cols; i++ ) {
            d_M[i] += alpha * M2[i];
        }
    }
}
void DenseSerialMatrix::setScalar( double alpha )
{
    for ( size_t i = 0; i < d_rows * d_cols; i++ )
        d_M[i]     = alpha;
}
void DenseSerialMatrix::zero() { memset( d_M, 0, d_rows * d_cols * sizeof( double ) ); }


/********************************************************
* Get/Set values                                        *
********************************************************/
void DenseSerialMatrix::addValuesByGlobalID(
    int num_rows, int num_cols, int *rows, int *cols, double *values )
{
    for ( int i = 0; i < num_rows; i++ ) {
        for ( int j = 0; j < num_cols; j++ ) {
            d_M[rows[i] + cols[j] * d_rows] += values[num_cols * i + j];
        }
    }
}
void DenseSerialMatrix::setValuesByGlobalID(
    int num_rows, int num_cols, int *rows, int *cols, double *values )
{
    for ( int i = 0; i < num_rows; i++ ) {
        for ( int j = 0; j < num_cols; j++ ) {
            d_M[rows[i] + cols[j] * d_rows] = values[num_cols * i + j];
        }
    }
}
void DenseSerialMatrix::addValueByGlobalID( int row, int col, double value )
{
    AMP_ASSERT( row < (int) d_rows && col < (int) d_cols );
    d_M[row + col * d_rows] += value;
}
void DenseSerialMatrix::setValueByGlobalID( int row, int col, double value )
{
    AMP_ASSERT( row < (int) d_rows && col < (int) d_cols );
    d_M[row + col * d_rows] = value;
}


/********************************************************
* Get values/rows by global id                          *
********************************************************/
void DenseSerialMatrix::getRowByGlobalID( int row,
                                          std::vector<unsigned int> &cols,
                                          std::vector<double> &values ) const
{
    AMP_ASSERT( row < (int) d_rows );
    cols.resize( d_cols );
    values.resize( d_cols );
    for ( size_t i = 0; i < d_cols; i++ ) {
        cols[i]   = i;
        values[i] = d_M[row + i * d_rows];
    }
}
double DenseSerialMatrix::getValueByGlobalID( int row, int col ) const
{
    return d_M[row + col * d_rows];
}
void DenseSerialMatrix::getValuesByGlobalID(
    int num_rows, int num_cols, int *rows, int *cols, double *values ) const
{
    for ( int i = 0; i < num_rows; i++ )
        for ( int j                  = 0; j < num_cols; j++ )
            values[i * num_cols + j] = d_M[rows[i] + cols[j] * d_rows];
}


/********************************************************
* Get/Set the diagonal                                  *
********************************************************/
Vector::shared_ptr DenseSerialMatrix::extractDiagonal( Vector::shared_ptr buf ) const
{
    AMP_ASSERT( d_cols == d_rows );
    Vector::shared_ptr out = buf;
    if ( buf == nullptr )
        out = this->getRightVector();
    AMP_ASSERT( out->getGlobalSize() == d_cols );
    double *y = new double[d_cols];
    for ( size_t i = 0; i < d_cols; i++ )
        y[i]       = d_M[i + i * d_rows];
    size_t *k      = new size_t[d_cols];
    for ( size_t i = 0; i < d_cols; i++ )
        k[i]       = i;
    out->setValuesByGlobalID( d_cols, k, y );
    delete[] y;
    delete[] k;
    return out;
}
void DenseSerialMatrix::setDiagonal( Vector::const_shared_ptr in )
{
    AMP_ASSERT( d_cols == d_rows );
    AMP_ASSERT( in->getGlobalSize() == d_rows );
    size_t *k = new size_t[d_rows];
    for ( size_t i = 0; i < d_rows; i++ )
        k[i]       = i;
    double *x      = new double[d_rows];
    in->getValuesByGlobalID( d_rows, k, x );
    for ( size_t i          = 0; i < d_rows; i++ )
        d_M[i + i * d_rows] = x[i];
    delete[] x;
}
void DenseSerialMatrix::setIdentity()
{
    AMP_ASSERT( d_cols == d_rows );
    memset( d_M, 0, d_rows * d_cols * sizeof( double ) );
    for ( size_t i          = 0; i < d_rows; i++ )
        d_M[i + i * d_rows] = 1.0;
}


/********************************************************
* Get the left/right vectors and DOFManagers            *
********************************************************/
Vector::shared_ptr DenseSerialMatrix::getRightVector() const
{
    return createVector( getRightDOFManager(), d_VariableRight );
}
Vector::shared_ptr DenseSerialMatrix::getLeftVector() const
{
    return createVector( getLeftDOFManager(), d_VariableLeft );
}
Discretization::DOFManager::shared_ptr DenseSerialMatrix::getRightDOFManager() const
{
    return d_DOFManagerRight;
}
Discretization::DOFManager::shared_ptr DenseSerialMatrix::getLeftDOFManager() const
{
    return d_DOFManagerLeft;
}


/********************************************************
* Compute the maximum column sum                        *
********************************************************/
double DenseSerialMatrix::L1Norm() const
{
    double norm = 0.0;
    for ( size_t j = 0; j < d_cols; j++ ) {
        double sum = 0.0;
        for ( size_t i = 0; i < d_rows; i++ )
            sum += fabs( d_M[i + j * d_rows] );
        norm = std::max( norm, sum );
    }
    return norm;
}

/********************************************************
* Multiply two matricies                                *
* result = this * other_op                              *
* C(N,M) = A(N,K)*B(K,M)
********************************************************/
void DenseSerialMatrix::multiply( Matrix::shared_ptr other_op, Matrix::shared_ptr &result )
{
    if ( this->numGlobalColumns() != other_op->numGlobalRows() )
        AMP_ERROR( "Inner matrix dimensions must agree" );
    size_t N = this->numGlobalRows();
    size_t K = this->numGlobalColumns();
    size_t M = other_op->numGlobalColumns();
    // Create the matrix parameters
    AMP::shared_ptr<AMP::LinearAlgebra::MatrixParameters> params(
        new AMP::LinearAlgebra::MatrixParameters(
            other_op->getRightDOFManager(), d_DOFManagerLeft, d_comm ) );
    params->d_VariableLeft  = d_VariableLeft;
    params->d_VariableRight = d_VariableRight;
    // Create the matrix
    AMP::shared_ptr<AMP::LinearAlgebra::DenseSerialMatrix> newMatrix(
        new AMP::LinearAlgebra::DenseSerialMatrix( params ) );
    result = newMatrix;
    memset( newMatrix->d_M, 0, N * M * sizeof( double ) );
    // Perform the muliplication
    if ( AMP::dynamic_pointer_cast<DenseSerialMatrix>( other_op ) == nullptr ) {
        // X is an unknown matrix type
        AMP_ERROR( "Not programmed yet" );
    } else {
        // We are dealing with all DenseSerialMatrix classes
        const double *A = d_M;
        const double *B = AMP::dynamic_pointer_cast<DenseSerialMatrix>( other_op )->d_M;
        double *C       = newMatrix->d_M;
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


} // LinearAlgebra namespace
} // AMP namespace
