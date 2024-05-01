#ifndef included_AMP_CSRMatrixData_hpp
#define included_AMP_CSRMatrixData_hpp

#include "AMP/AMP_TPLs.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/CSRMatrixParameters.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Utilities.h"

#ifdef AMP_USE_UMPIRE
    #include "umpire/Allocator.hpp"
    #include "umpire/ResourceManager.hpp"
#endif

#include <memory>
#include <numeric>
#include <type_traits>

namespace AMP::LinearAlgebra {


/********************************************************
 * Constructors/Destructor                               *
 ********************************************************/
template<typename Policy>
CSRMatrixData<Policy>::CSRMatrixData()
{
    AMPManager::incrementResource( "CSRMatrixData" );
}

template<typename Policy>
static auto NNZ( typename Policy::lidx_t N, typename Policy::lidx_t *nnz_per_row ) ->
    typename Policy::lidx_t
{
    AMP_ASSERT( AMP::Utilities::getMemoryType( nnz_per_row ) < AMP::Utilities::MemoryType::device );
    return std::accumulate( nnz_per_row, nnz_per_row + N, 0 );
}

template<typename Policy>
CSRMatrixData<Policy>::CSRMatrixData( std::shared_ptr<MatrixParametersBase> params )
    : MatrixData( params )
{
    AMPManager::incrementResource( "CSRMatrixData" );
    auto csrParams = std::dynamic_pointer_cast<CSRMatrixParameters<Policy>>( d_pParameters );
    auto matParams = std ::dynamic_pointer_cast<MatrixParameters>( d_pParameters );

    d_memory_location = d_pParameters->d_memory_location;

    if ( csrParams ) {

        AMP_ERROR( "Should this even be allowed, need to know more about DOFManagers. " );

    } else if ( matParams ) {

        // for now all matrix parameter data is assumed to be on host

        d_leftDOFManager  = matParams->getLeftDOFManager();
        d_rightDOFManager = matParams->getRightDOFManager();
        AMP_ASSERT( d_leftDOFManager && d_rightDOFManager );

        d_is_square = ( d_leftDOFManager->numGlobalDOF() == d_rightDOFManager->numGlobalDOF() );
        d_first_row = d_leftDOFManager->beginDOF();
        d_last_row  = d_leftDOFManager->endDOF();
        d_first_col = d_rightDOFManager->beginDOF();
        d_last_col  = d_rightDOFManager->endDOF();

	d_diagMatrix = std::make_shared<CSRSerialMatrixData<Policy>>( params, true );
	d_offDiagMatrix = std::make_shared<CSRSerialMatrixData<Policy>>( params, false );

	d_nnz = d_diagMatrix->numberOfNonZeros() + d_offDiagMatrix->numberOfNonZeros();

    } else {
        AMP_ERROR( "Check supplied MatrixParameter object" );
    }
}

template<typename Policy>
CSRMatrixData<Policy>::~CSRMatrixData()
{
    AMPManager::decrementResource( "CSRMatrixData" );
}

template<typename Policy>
std::shared_ptr<MatrixData> CSRMatrixData<Policy>::cloneMatrixData() const
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
std::shared_ptr<MatrixData> CSRMatrixData<Policy>::transpose() const
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixData<Policy>::extractDiagonal( std::shared_ptr<Vector> buf ) const
{
    // Forward to diagonal matrix
    d_diagMatrix->extractDiagonal( buf );
}
template<typename Policy>
void CSRMatrixData<Policy>::getRowByGlobalID( size_t row,
                                              std::vector<size_t> &cols,
                                              std::vector<double> &vals ) const
{
    // Get portion of row from diagonal matrix
    d_diagMatrix->getRowByGlobalID( row, cols, vals );
    // Get portion from off diagonal and append
    std::vector<size_t> od_cols;
    std::vector<double> od_vals;
    d_offDiagMatrix->getRowByGlobalID( row, od_cols, od_vals );
    cols.insert( cols.end(), od_cols.begin(), od_cols.end() );
    vals.insert( vals.end(), od_vals.begin(), od_vals.end() );
}

template<typename Policy>
void CSRMatrixData<Policy>::addValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, void *vals, const typeID &id )
{
    // Internal sparse matrices will ignore entries they don't own
    // this can just be forwarded
    d_diagMatrix->addValuesByGlobalID( num_rows, num_cols, rows, cols, vals, id);
    d_offDiagMatrix->addValuesByGlobalID( num_rows, num_cols, rows, cols, vals, id);
}

template<typename Policy>
void CSRMatrixData<Policy>::setValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, void *vals, const typeID &id )
{
    // Internal sparse matrices will ignore entries they don't own
    // this can just be forwarded
    d_diagMatrix->setValuesByGlobalID( num_rows, num_cols, rows, cols, vals, id);
    d_offDiagMatrix->setValuesByGlobalID( num_rows, num_cols, rows, cols, vals, id);
}

template<typename Policy>
void CSRMatrixData<Policy>::getValuesByGlobalID( size_t num_rows,
                                                 size_t num_cols,
                                                 size_t *rows,
                                                 size_t *cols,
                                                 void *values,
                                                 const typeID &id ) const
{
    // Forward to internal matrices, nothing will happen if not found
    d_diagMatrix->getValuesByGlobalID( num_rows, num_cols, rows, cols, values, id );
    d_offDiagMatrix->getValuesByGlobalID( num_rows, num_cols, rows, cols, values, id );
}

template<typename Policy>
void CSRMatrixData<Policy>::setOtherData( std::map<gidx_t, std::map<gidx_t, scalar_t>> &other_data,
                                          AMP::LinearAlgebra::ScatterType t )
{
    AMP_MPI comm   = getComm();
    auto ndxLen    = other_data.size();
    auto totNdxLen = comm.sumReduce( ndxLen );
    if ( totNdxLen == 0 ) {
        return;
    }
    auto dataLen = 0;
    auto cur_row = other_data.begin();
    while ( cur_row != other_data.end() ) {
        dataLen += cur_row->second.size();
        ++cur_row;
    }
    std::vector<gidx_t> rows( dataLen + 1 );   // Add one to have the new work
    std::vector<gidx_t> cols( dataLen + 1 );   // Add one to have the new work
    std::vector<scalar_t> data( dataLen + 1 ); // Add one to have the new work
    size_t cur_ptr = 0;
    cur_row        = other_data.begin();
    while ( cur_row != other_data.end() ) {
        auto cur_elem = cur_row->second.begin();
        while ( cur_elem != cur_row->second.end() ) {
            rows[cur_ptr] = cur_row->first;
            cols[cur_ptr] = cur_elem->first;
            data[cur_ptr] = cur_elem->second;
            ++cur_ptr;
            ++cur_elem;
        }
        ++cur_row;
    }

    auto totDataLen = comm.sumReduce( dataLen );

    std::vector<gidx_t> aggregateRows( totDataLen );
    std::vector<gidx_t> aggregateCols( totDataLen );
    std::vector<scalar_t> aggregateData( totDataLen );

    comm.allGather( rows.data(), dataLen, aggregateRows.data() );
    comm.allGather( cols.data(), dataLen, aggregateCols.data() );
    comm.allGather( data.data(), dataLen, aggregateData.data() );

    if ( t == AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD ) {
        for ( int i = 0; i != totDataLen; i++ ) {
            if ( ( aggregateRows[i] >= d_first_row ) && ( aggregateRows[i] < d_last_row ) ) {
                addValuesByGlobalID( 1u,
                                     1u,
                                     (size_t *) &aggregateRows[i],
                                     (size_t *) &aggregateCols[i],
                                     &aggregateData[i],
                                     getTypeID<scalar_t>() );
            }
        }
    } else {

        if ( t == AMP::LinearAlgebra::ScatterType::CONSISTENT_SET ) {
            for ( int i = 0; i != totDataLen; i++ ) {
                if ( ( aggregateRows[i] >= d_first_row ) && ( aggregateRows[i] < d_last_row ) ) {
                    setValuesByGlobalID( 1u,
                                         1u,
                                         (size_t *) &aggregateRows[i],
                                         (size_t *) &aggregateCols[i],
                                         &aggregateData[i],
                                         getTypeID<scalar_t>() );
                }
            }
        }
    }

    other_data.clear();
}

template<typename Policy>
void CSRMatrixData<Policy>::makeConsistent( AMP::LinearAlgebra::ScatterType t )
{
    if ( t == AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD )
        setOtherData( d_other_data, AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
    else
        setOtherData( d_ghost_data, AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}


template<typename Policy>
std::vector<size_t> CSRMatrixData<Policy>::getColumnIDs( size_t row ) const
{
    AMP_INSIST( d_diagMatrix, "diag matrix must exist" );
    std::vector<size_t> cols = d_diagMatrix->getColumnIDs( row );
    if ( hasOffDiag() ) {
        std::vector<size_t> od_cols = d_offDiagMatrix->getColumnIDs( row );
	cols.insert( cols.end(), od_cols.begin(), od_cols.end() );
    }
    return cols;
}

template<typename Policy>
std::shared_ptr<Discretization::DOFManager> CSRMatrixData<Policy>::getRightDOFManager() const
{
    return d_rightDOFManager;
}

template<typename Policy>
std::shared_ptr<Discretization::DOFManager> CSRMatrixData<Policy>::getLeftDOFManager() const
{
    return d_leftDOFManager;
}

/********************************************************
 * Get the number of rows/columns in the matrix          *
 ********************************************************/
template<typename Policy>
size_t CSRMatrixData<Policy>::numLocalRows() const
{
    return static_cast<size_t>( d_last_row - d_first_row );
}

template<typename Policy>
size_t CSRMatrixData<Policy>::numGlobalRows() const
{
    AMP_ASSERT( d_leftDOFManager );
    return d_leftDOFManager->numGlobalDOF();
}

template<typename Policy>
size_t CSRMatrixData<Policy>::numLocalColumns() const
{
    return static_cast<size_t>( d_last_col - d_first_col );
}

template<typename Policy>
size_t CSRMatrixData<Policy>::numGlobalColumns() const
{
    AMP_ASSERT( d_rightDOFManager );
    return d_rightDOFManager->numGlobalDOF();
}

/********************************************************
 * Get iterators                                         *
 ********************************************************/
template<typename Policy>
size_t CSRMatrixData<Policy>::beginRow() const
{
    return static_cast<size_t>( d_first_row );
}

template<typename Policy>
size_t CSRMatrixData<Policy>::endRow() const
{
    return static_cast<size_t>( d_last_row );
}

} // namespace AMP::LinearAlgebra

#endif
