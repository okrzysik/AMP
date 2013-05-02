#include "matrices/DenseSerialMatrix.h"
#include <stdio.h>
#include <string.h>


namespace AMP {
namespace LinearAlgebra {


/********************************************************
* Constructor/Destructor                                *
********************************************************/
DenseSerialMatrix::DenseSerialMatrix( MatrixParameters::shared_ptr  params ):
    Matrix(params)
{
    d_VariableLeft = params->d_VariableLeft;
    d_VariableRight = params->d_VariableRight;
    d_DOFManagerLeft = params->getLeftDOFManager();
    d_DOFManagerRight = params->getRightDOFManager();
    d_rows = params->getGlobalNumberOfRows();
    d_cols = params->getGlobalNumberOfColumns();
    d_M = new double[d_rows*d_cols];
    memset(d_M,0,d_rows*d_cols*sizeof(double));
}
DenseSerialMatrix::~DenseSerialMatrix()
{
    delete [] d_M;
}


/********************************************************
* Copy/transpose the matrix                             *
********************************************************/
Matrix::shared_ptr  DenseSerialMatrix::transpose() const
{
    AMP_ERROR("Not implimented yet");
    return Matrix::shared_ptr();
}
Matrix::shared_ptr  DenseSerialMatrix::cloneMatrix() const
{
    AMP_ERROR("Not implimented yet");
    return Matrix::shared_ptr();
}


/********************************************************
* Matrix-vector multiplication                          *
********************************************************/
void DenseSerialMatrix::mult( 
    AMP::LinearAlgebra::Vector::const_shared_ptr in, 
    AMP::LinearAlgebra::Vector::shared_ptr out )
{
    AMP_ASSERT(in->getGlobalSize()==d_cols);
    AMP_ASSERT(out->getGlobalSize()==d_rows);
    size_t *k = new size_t[std::max(d_cols,d_rows)];
    for (size_t i=0; i<std::max(d_cols,d_rows); i++)
        k[i] = i;
    // Get x
    double *x = new double[d_cols];
    in->getValuesByGlobalID(d_cols,k,x);    
    // Initialize y
    double *y = new double[d_rows];
    memset(y,0,d_rows*sizeof(double));
    // Perform y = M*x
    for (size_t j=0; j<d_cols; j++) {
        for (size_t i=0; i<d_rows; i++)
            y[i] += d_M[i+j*d_rows]*x[j];
    }
    // Save y
    out->setValuesByGlobalID(d_rows,k,y);    
    delete [] x;
    delete [] y;
    delete [] k;
}
void DenseSerialMatrix::multTranspose( 
    AMP::LinearAlgebra::Vector::const_shared_ptr in, 
    AMP::LinearAlgebra::Vector::shared_ptr out )
{
    AMP_ASSERT(in->getGlobalSize()==d_rows);
    AMP_ASSERT(out->getGlobalSize()==d_cols);
    size_t *k = new size_t[std::max(d_cols,d_rows)];
    for (size_t i=0; i<std::max(d_cols,d_rows); i++)
        k[i] = i;
    // Get x
    double *x = new double[d_rows];
    in->getValuesByGlobalID(d_rows,k,x);    
    // Initialize y
    double *y = new double[d_cols];
    memset(y,0,d_cols*sizeof(double));
    // Perform y = M*x
    for (size_t j=0; j<d_cols; j++) {
        for (size_t i=0; i<d_rows; i++)
            y[j] += d_M[i+j*d_rows]*x[i];
    }
    // Save y
    out->setValuesByGlobalID(d_cols,k,y);    
    delete [] x;
    delete [] y;
    delete [] k;
}


/********************************************************
* Scale/axpy/setScalar                                  *
********************************************************/
void DenseSerialMatrix::scale( double alpha )
{
    for (size_t i=0; i<d_rows*d_cols; i++)
        d_M[i] *= alpha;
}
void DenseSerialMatrix::axpy( double alpha, const Matrix &X )
{
    AMP_ERROR("Not implimented yet");
}
void DenseSerialMatrix::setScalar( double alpha )
{
    for (size_t i=0; i<d_rows*d_cols; i++)
        d_M[i] = alpha;
}


/********************************************************
* Get/Set values                                        *
********************************************************/
void DenseSerialMatrix::addValuesByGlobalID( int   num_rows ,
                                          int   num_cols ,
                                          int  *rows ,
                                          int  *cols ,
                                          double  *values )
{
    for (int i=0; i<num_rows; i++) {
        for (int j=0; j<num_cols; j++) {
            d_M[rows[i]*cols[j]*d_rows] += values[num_cols*i+j];
        }
    }
}
void DenseSerialMatrix::setValuesByGlobalID( int   num_rows ,
                                          int   num_cols ,
                                          int  *rows ,
                                          int  *cols ,
                                          double  *values )
{
    for (int i=0; i<num_rows; i++) {
        for (int j=0; j<num_cols; j++) {
            d_M[rows[i]*cols[j]*d_rows] = values[num_cols*i+j];
        }
    }
}
void DenseSerialMatrix::addValueByGlobalID( int row , int col , double value )
{
    AMP_ASSERT(row<(int)d_rows&&col<(int)d_cols);
    d_M[row+col*d_rows] += value;
}
void DenseSerialMatrix::setValueByGlobalID( int row , int col , double value )
{
    AMP_ASSERT(row<(int)d_rows&&col<(int)d_cols);
    d_M[row+col*d_rows] = value;
}


/********************************************************
* Get a row                                             *
********************************************************/
void DenseSerialMatrix::getRowByGlobalID( int row, std::vector<unsigned int> &cols , 
                                                  std::vector<double> &values ) const
{
    AMP_ASSERT(row<(int)d_rows);
    cols.resize(d_cols);
    values.resize(d_cols);
    for (size_t i=0; i<d_cols; i++) {
        cols[i] = i;
        values[i] = d_M[row+i*d_rows];
    }
}


/********************************************************
* Get/Set the diagonal                                  *
********************************************************/
Vector::shared_ptr  DenseSerialMatrix::extractDiagonal( Vector::shared_ptr buf )
{
    AMP_ASSERT(d_cols==d_rows);
    Vector::shared_ptr out = buf;
    if ( buf == NULL )
        out = this->getRightVector();
    AMP_ASSERT(out->getGlobalSize()==d_cols);
    double *y = new double[d_cols];
    for (size_t i=0; i<d_cols; i++)
        y[i] = d_M[i+i*d_rows];
    size_t *k = new size_t[d_cols];
    for (size_t i=0; i<d_cols; i++)
        k[i] = i;
    out->setValuesByGlobalID(d_cols,k,y);
    delete [] y;
    delete [] k;
    return out;
}
void DenseSerialMatrix::setDiagonal( const Vector::shared_ptr &in )
{
    AMP_ASSERT(d_cols==d_rows);
    AMP_ASSERT(in->getGlobalSize()==d_rows);
    size_t *k = new size_t[d_rows];
    for (size_t i=0; i<d_rows; i++)
        k[i] = i;
    double *x = new double[d_rows];
    in->getValuesByGlobalID(d_rows,k,x);
    for (size_t i=0; i<d_rows; i++)
        d_M[i+i*d_rows] = x[i];
    delete [] x;
}


/********************************************************
* Get the left/right vectors and DOFManagers            *
********************************************************/
Vector::shared_ptr DenseSerialMatrix::getRightVector()
{
    AMP_ERROR("Not implimented yet");
    return Vector::shared_ptr();
}
Vector::shared_ptr DenseSerialMatrix::getLeftVector()
{
    AMP_ERROR("Not implimented yet");
    return Vector::shared_ptr();
}
Discretization::DOFManager::shared_ptr DenseSerialMatrix::getRightDOFManager()
{
    return d_DOFManagerRight;
}
Discretization::DOFManager::shared_ptr DenseSerialMatrix::getLeftDOFManager()
{
    return d_DOFManagerLeft;
}


/********************************************************
* Compute the maximum column sum                        *
********************************************************/
double DenseSerialMatrix::L1Norm() const
{
    AMP_ERROR("Not implimented yet");
    return 0.0;
}

/********************************************************
* Multiply two matricies                                *
********************************************************/
void DenseSerialMatrix::multiply( Matrix::shared_ptr other_op, Matrix::shared_ptr &result )
{
    AMP_ERROR("Not implimented yet");
}



} // LinearAlgebra namespace
} // AMP namespace

