#include "AMP/matrices/data/DenseSerialMatrixData.h"
#include "AMP/vectors/VectorBuilder.h"
#include <cstdio>
#include <cstring>

#include <numeric>

namespace AMP::LinearAlgebra {


/********************************************************
 * Constructor/Destructor                                *
 ********************************************************/
DenseSerialMatrixData::DenseSerialMatrixData( std::shared_ptr<MatrixParameters> params )
    : MatrixData( params ),
      d_VariableLeft( params->d_VariableLeft ),
      d_VariableRight( params->d_VariableRight ),
      d_DOFManagerLeft( params->getLeftDOFManager() ),
      d_DOFManagerRight( params->getRightDOFManager() ),
      d_rows( params->getGlobalNumberOfRows() ),
      d_cols( params->getGlobalNumberOfColumns() )
{
    d_M = new double[d_rows * d_cols];
    memset( d_M, 0, d_rows * d_cols * sizeof( double ) );
}

DenseSerialMatrixData::~DenseSerialMatrixData() { delete[] d_M; }

std::shared_ptr<MatrixData> DenseSerialMatrixData::cloneMatrixData() const
{
    // Create the matrix parameters
    auto params = std::make_shared<AMP::LinearAlgebra::MatrixParameters>(
        d_DOFManagerLeft, d_DOFManagerRight, getComm() );
    params->d_VariableLeft  = d_VariableLeft;
    params->d_VariableRight = d_VariableRight;
    // Create the matrix
    auto newMatrixData = std::make_shared<AMP::LinearAlgebra::DenseSerialMatrixData>( params );
    double *M2         = newMatrixData->d_M;
    memcpy( M2, d_M, d_cols * d_rows * sizeof( double ) );
    return newMatrixData;
}

/********************************************************
 * Get/Set values                                        *
 ********************************************************/
void DenseSerialMatrixData::addValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values )
{
    if ( num_rows == 1 && num_cols == 1 ) {
        d_M[rows[0] + cols[0] * d_rows] += values[0];
    } else {
        for ( size_t i = 0; i < num_rows; i++ ) {
            for ( size_t j = 0; j < num_cols; j++ ) {
                d_M[rows[i] + cols[j] * d_rows] += values[num_cols * i + j];
            }
        }
    }
}
void DenseSerialMatrixData::setValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values )
{
    if ( num_rows == 1 && num_cols == 1 ) {
        d_M[rows[0] + cols[0] * d_rows] = values[0];
    } else {
        for ( size_t i = 0; i < num_rows; i++ ) {
            for ( size_t j = 0; j < num_cols; j++ ) {
                d_M[rows[i] + cols[j] * d_rows] = values[num_cols * i + j];
            }
        }
    }
}
void DenseSerialMatrixData::getValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) const
{
    if ( num_rows == 1 && num_cols == 1 ) {
        values[0] = d_M[rows[0] + cols[0] * d_rows];
    } else {
        for ( size_t i = 0; i < num_rows; i++ )
            for ( size_t j = 0; j < num_cols; j++ )
                values[i * num_cols + j] = d_M[rows[i] + cols[j] * d_rows];
    }
}

/********************************************************
 * Get values/rows by global id                          *
 ********************************************************/
void DenseSerialMatrixData::getRowByGlobalID( size_t row,
                                              std::vector<size_t> &cols,
                                              std::vector<double> &values ) const
{
    AMP_ASSERT( row < d_rows );
    cols.resize( d_cols );
    values.resize( d_cols );
    for ( size_t i = 0; i < d_cols; i++ ) {
        cols[i]   = i;
        values[i] = d_M[row + i * d_rows];
    }
}

/********************************************************
 * Get column indices by global id                       *
 ********************************************************/
std::vector<size_t> DenseSerialMatrixData::getColumnIDs( size_t row ) const
{
    AMP_ASSERT( row < d_rows );

    std::vector<size_t> cols( d_cols );
    std::iota( cols.begin(), cols.end(), 0 );

    return cols;
}

/********************************************************
 * Get the left/right DOFManagers                       *
 ********************************************************/

Discretization::DOFManager::shared_ptr DenseSerialMatrixData::getRightDOFManager() const
{
    return d_DOFManagerRight;
}
Discretization::DOFManager::shared_ptr DenseSerialMatrixData::getLeftDOFManager() const
{
    return d_DOFManagerLeft;
}

} // namespace AMP::LinearAlgebra
