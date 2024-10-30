#ifndef included_AMP_CSRMatrixData_hpp
#define included_AMP_CSRMatrixData_hpp

#include "AMP/AMP_TPLs.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/CSRMatrixParameters.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Utilities.h"

#include <type_traits>

namespace AMP::LinearAlgebra {

/********************************************************
 * Constructors/Destructor                              *
 ********************************************************/
template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::CSRMatrixData()
    : d_memory_location( AMP::Utilities::getAllocatorMemoryType<Allocator>() )
{
    AMPManager::incrementResource( "CSRMatrixData" );
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::CSRMatrixData(
    std::shared_ptr<MatrixParametersBase> params )
    : MatrixData( params ), d_memory_location( AMP::Utilities::getAllocatorMemoryType<Allocator>() )
{
    AMPManager::incrementResource( "CSRMatrixData" );
    auto csrParams = std::dynamic_pointer_cast<CSRMatrixParameters<Policy>>( d_pParameters );
    auto matParams = std ::dynamic_pointer_cast<MatrixParameters>( d_pParameters );

    AMP_INSIST(
        d_memory_location != AMP::Utilities::MemoryType::device,
        "CSRMatrixData and CSRSerialMatrixData do not support pure-device memory locations yet" );

    if ( csrParams ) {

        d_first_row = csrParams->d_first_row;
        d_last_row  = csrParams->d_last_row;
        d_first_col = csrParams->d_first_col;
        d_last_col  = csrParams->d_last_col;

        // Construct on/off diag blocks
        d_diag_matrix = std::make_shared<DiagMatrixData>(
            params, d_memory_location, d_first_row, d_last_row, d_first_col, d_last_col, true );
        d_offd_matrix = std::make_shared<OffdMatrixData>(
            params, d_memory_location, d_first_row, d_last_row, d_first_col, d_last_col, false );

        // Make DOF managers
        std::vector<size_t> remoteDOFsRight;
        d_offd_matrix->getColumnMap( remoteDOFsRight );
        d_rightDOFManager = std::make_shared<Discretization::DOFManager>(
            d_last_row - d_first_row, getComm(), remoteDOFsRight );
        // Finding ghosts for leftDM hard to do and not currently needed
        // should think about how to approach it though
        d_leftDOFManager =
            std::make_shared<Discretization::DOFManager>( d_last_col - d_first_col, getComm() );

    } else if ( matParams ) {

        d_leftDOFManager  = matParams->getLeftDOFManager();
        d_rightDOFManager = matParams->getRightDOFManager();
        AMP_ASSERT( d_leftDOFManager && d_rightDOFManager );
        d_first_row = d_leftDOFManager->beginDOF();
        d_last_row  = d_leftDOFManager->endDOF();
        d_first_col = d_rightDOFManager->beginDOF();
        d_last_col  = d_rightDOFManager->endDOF();

        // send params forward to the on/off diagonal blocks
        d_diag_matrix = std::make_shared<DiagMatrixData>(
            params, d_memory_location, d_first_row, d_last_row, d_first_col, d_last_col, true );
        d_offd_matrix = std::make_shared<OffdMatrixData>(
            params, d_memory_location, d_first_row, d_last_row, d_first_col, d_last_col, false );
        d_nnz = d_diag_matrix->d_nnz + d_offd_matrix->d_nnz;

    } else {
        AMP_ERROR( "Check supplied MatrixParameters object" );
    }

    // get total nnz count
    d_nnz = d_diag_matrix->d_nnz + d_offd_matrix->d_nnz;

    d_is_square = ( d_leftDOFManager->numGlobalDOF() == d_rightDOFManager->numGlobalDOF() );
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::~CSRMatrixData()
{
    AMPManager::decrementResource( "CSRMatrixData" );
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
std::shared_ptr<MatrixData>
CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::cloneMatrixData() const
{
    std::shared_ptr<CSRMatrixData> cloneData;

    cloneData =
        std::make_shared<CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>>();

    cloneData->d_is_square       = d_is_square;
    cloneData->d_first_row       = d_first_row;
    cloneData->d_last_row        = d_last_row;
    cloneData->d_first_col       = d_first_col;
    cloneData->d_last_col        = d_last_col;
    cloneData->d_nnz             = d_nnz;
    cloneData->d_leftDOFManager  = d_leftDOFManager;
    cloneData->d_rightDOFManager = d_rightDOFManager;
    cloneData->d_pParameters     = d_pParameters;

    cloneData->d_diag_matrix = d_diag_matrix->cloneMatrixData();
    cloneData->d_offd_matrix = d_offd_matrix->cloneMatrixData();

    return cloneData;
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
std::shared_ptr<MatrixData>
CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::transpose() const
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::getRowByGlobalID(
    size_t row, std::vector<size_t> &cols, std::vector<double> &vals ) const
{
    AMP_DEBUG_INSIST( row >= static_cast<size_t>( d_first_row ) &&
                          row < static_cast<size_t>( d_last_row ),
                      "row must be owned by rank" );

    AMP_DEBUG_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                      "CSRMatrixData::getRowByGlobalID not implemented for device memory" );

    auto local_row = row - d_first_row;

    // Get portion of row from diagonal matrix
    d_diag_matrix->getRowByGlobalID( local_row, cols, vals );

    // Get portion from off diagonal and append
    std::vector<size_t> od_cols;
    std::vector<double> od_vals;
    d_offd_matrix->getRowByGlobalID( local_row, od_cols, od_vals );
    cols.insert( cols.end(), od_cols.begin(), od_cols.end() );
    vals.insert( vals.end(), od_vals.begin(), od_vals.end() );
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::getValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, void *values, const typeID &id )
    const
{
    AMP_DEBUG_INSIST( getTypeID<scalar_t>() == id,
                      "CSRMatrixData::getValuesByGlobalID called with inconsistent typeID" );

    AMP_DEBUG_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                      "CSRMatrixData::getValuesByGlobalID not implemented for device memory" );

    if ( num_rows == 1 && num_cols == 1 ) {

        const auto local_row = rows[0] - d_first_row;
        // Forward to internal matrices, nothing will happen if not found
        d_diag_matrix->getValuesByGlobalID( local_row, cols[0], values, id );
        d_offd_matrix->getValuesByGlobalID( local_row, cols[0], values, id );
    } else {
        AMP_ERROR( "CSRSerialMatrixData::getValuesByGlobalID not implemented for num_rows > 1 || "
                   "num_cols > 1" );
    }
}

// The two getValues functions above can directly forward to the diag and off diag blocks
// The addValuesByGlobalID and setValuesByGlobalID functions can't do this since
// they need to also handle the other_data case
template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::addValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, void *vals, const typeID &id )
{
    AMP_DEBUG_INSIST( getTypeID<scalar_t>() == id,
                      "CSRMatrixData::addValuesByGlobalID called with inconsistent typeID" );

    AMP_DEBUG_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                      "CSRMatrixData::addValuesByGlobalID not implemented for device memory" );

    auto values = reinterpret_cast<const scalar_t *>( vals );

    for ( size_t i = 0u; i != num_rows; i++ ) {
        if ( rows[i] >= static_cast<size_t>( d_first_row ) &&
             rows[i] < static_cast<size_t>( d_last_row ) ) {

            // Forward single row to diag and off diag blocks
            // auto lcols = &cols[num_cols * i];
            const auto local_row = rows[i] - d_first_row;
            auto lvals           = &values[num_cols * i];
            d_diag_matrix->addValuesByGlobalID( num_cols, local_row, cols, lvals, id );
            d_offd_matrix->addValuesByGlobalID( num_cols, local_row, cols, lvals, id );
        } else {
            for ( size_t icol = 0; icol < num_cols; ++icol ) {
                d_other_data[rows[i]][cols[icol]] += values[num_cols * i + icol];
            }
        }
    }
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::setValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, void *vals, const typeID &id )
{
    AMP_DEBUG_INSIST( getTypeID<scalar_t>() == id,
                      "CSRMatrixData::setValuesByGlobalID called with inconsistent typeID" );

    AMP_DEBUG_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                      "CSRMatrixData::setValuesByGlobalID not implemented for device memory" );

    auto values = reinterpret_cast<const scalar_t *>( vals );

    for ( size_t i = 0u; i != num_rows; i++ ) {

        if ( rows[i] >= static_cast<size_t>( d_first_row ) &&
             rows[i] < static_cast<size_t>( d_last_row ) ) {

            // Forward single row to diag and off diag blocks
            // auto lcols = &cols[num_cols * i];
            const auto local_row = rows[i] - d_first_row;
            auto lvals           = &values[num_cols * i];
            d_diag_matrix->setValuesByGlobalID( num_cols, local_row, cols, lvals, id );
            d_offd_matrix->setValuesByGlobalID( num_cols, local_row, cols, lvals, id );
        } else {
            for ( size_t icol = 0; icol < num_cols; ++icol ) {
                d_ghost_data[rows[i]][cols[icol]] = values[num_cols * i + icol];
            }
        }
    }
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
std::vector<size_t>
CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::getColumnIDs( size_t row ) const
{
    AMP_DEBUG_INSIST( row >= static_cast<size_t>( d_first_row ) &&
                          row < static_cast<size_t>( d_last_row ),
                      "CSRMatrixData::getColumnIDs row must be owned by rank" );

    AMP_DEBUG_INSIST( d_diag_matrix, "CSRMatrixData::getColumnIDs diag matrix must exist" );

    AMP_DEBUG_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                      "CSRMatrixData::getColumnIDs not implemented for device memory" );

    auto local_row              = row - d_first_row;
    std::vector<size_t> cols    = d_diag_matrix->getColumnIDs( local_row );
    std::vector<size_t> od_cols = d_offd_matrix->getColumnIDs( local_row );
    cols.insert( cols.end(), od_cols.begin(), od_cols.end() );
    return cols;
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::setOtherData(
    std::map<gidx_t, std::map<gidx_t, scalar_t>> &other_data, AMP::LinearAlgebra::ScatterType t )
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

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
void CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::makeConsistent(
    AMP::LinearAlgebra::ScatterType t )
{
    if ( t == AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD )
        setOtherData( d_other_data, AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
    else
        setOtherData( d_ghost_data, AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
std::shared_ptr<Discretization::DOFManager>
CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::getRightDOFManager() const
{
    return d_rightDOFManager;
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
std::shared_ptr<Discretization::DOFManager>
CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::getLeftDOFManager() const
{
    return d_leftDOFManager;
}

/********************************************************
 * Get the number of rows/columns in the matrix          *
 ********************************************************/
template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
size_t CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::numLocalRows() const
{
    return static_cast<size_t>( d_last_row - d_first_row );
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
size_t CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::numGlobalRows() const
{
    AMP_ASSERT( d_leftDOFManager );
    return d_leftDOFManager->numGlobalDOF();
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
size_t CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::numLocalColumns() const
{
    return static_cast<size_t>( d_last_col - d_first_col );
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
size_t CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::numGlobalColumns() const
{
    AMP_ASSERT( d_rightDOFManager );
    return d_rightDOFManager->numGlobalDOF();
}

/********************************************************
 * Get iterators                                         *
 ********************************************************/
template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
size_t CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::beginRow() const
{
    return static_cast<size_t>( d_first_row );
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
size_t CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData>::endRow() const
{
    return static_cast<size_t>( d_last_row );
}

} // namespace AMP::LinearAlgebra

#endif
