#ifndef included_AMP_CSRMatrixData_hpp
#define included_AMP_CSRMatrixData_hpp

#include "AMP/AMP_TPLs.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/AMPCSRMatrixParameters.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/matrices/MatrixParametersBase.h"
#include "AMP/matrices/RawCSRMatrixParameters.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Utilities.h"

#include <type_traits>

namespace AMP::LinearAlgebra {

/********************************************************
 * Constructors/Destructor                              *
 ********************************************************/
template<typename Policy, class Allocator, class DiagMatrixData>
CSRMatrixData<Policy, Allocator, DiagMatrixData>::CSRMatrixData()
    : d_memory_location( AMP::Utilities::getAllocatorMemoryType<Allocator>() )
{
    AMPManager::incrementResource( "CSRMatrixData" );
}

template<typename Policy, class Allocator, class DiagMatrixData>
CSRMatrixData<Policy, Allocator, DiagMatrixData>::CSRMatrixData(
    std::shared_ptr<MatrixParametersBase> params )
    : MatrixData( params ), d_memory_location( AMP::Utilities::getAllocatorMemoryType<Allocator>() )
{
    AMPManager::incrementResource( "CSRMatrixData" );

    AMP_INSIST(
        d_memory_location != AMP::Utilities::MemoryType::device,
        "CSRMatrixData and CSRSerialMatrixData do not support pure-device memory locations yet" );

    // Figure out what kind of parameters object we have
    // Note: matParams always true if ampCSRParams is by inheritance
    auto rawCSRParams = std::dynamic_pointer_cast<RawCSRMatrixParameters<Policy>>( d_pParameters );
    auto ampCSRParams = std::dynamic_pointer_cast<AMPCSRMatrixParameters<Policy>>( d_pParameters );
    auto matParams    = std ::dynamic_pointer_cast<MatrixParameters>( d_pParameters );

    if ( rawCSRParams ) {

        // Simplest initialization, extract row/column bounds and pass through to diag/offd
        d_first_row = rawCSRParams->d_first_row;
        d_last_row  = rawCSRParams->d_last_row;
        d_first_col = rawCSRParams->d_first_col;
        d_last_col  = rawCSRParams->d_last_col;

        // Construct on/off diag blocks
        d_diag_matrix = std::make_shared<DiagMatrixData>(
            params, d_memory_location, d_first_row, d_last_row, d_first_col, d_last_col, true );
        d_offd_matrix = std::make_shared<DiagMatrixData>(
            params, d_memory_location, d_first_row, d_last_row, d_first_col, d_last_col, false );

        d_leftDOFManager  = nullptr;
        d_rightDOFManager = nullptr;
        d_leftCommList    = nullptr;
        d_rightCommList   = nullptr;

    } else if ( matParams ) {
        // Pull out DOFManagers and CommunicationLists, set row/col bounds
        d_leftDOFManager  = matParams->getLeftDOFManager();
        d_rightDOFManager = matParams->getRightDOFManager();
        AMP_ASSERT( d_leftDOFManager && d_rightDOFManager );
        d_first_row     = d_leftDOFManager->beginDOF();
        d_last_row      = d_leftDOFManager->endDOF();
        d_first_col     = d_rightDOFManager->beginDOF();
        d_last_col      = d_rightDOFManager->endDOF();
        d_leftCommList  = matParams->getLeftCommList();
        d_rightCommList = matParams->getRightCommList();

        // Construct on/off diag blocks
        d_diag_matrix = std::make_shared<DiagMatrixData>(
            params, d_memory_location, d_first_row, d_last_row, d_first_col, d_last_col, true );
        d_offd_matrix = std::make_shared<DiagMatrixData>(
            params, d_memory_location, d_first_row, d_last_row, d_first_col, d_last_col, false );

        // If, more specifically, have ampCSRParams then blocks are not yet
        // filled. This consolidates calls to getRow{NNZ,Cols} for both blocks
        if ( ampCSRParams && ampCSRParams->getRowNNZFunction() ) {
            AMP_INSIST( ampCSRParams->getRowColsFunction(),
                        "ampCSRParams->getRowColsFunction() must give valid function" );

            // number of local rows
            const lidx_t nrows = static_cast<lidx_t>( d_last_row - d_first_row );

            // pull out row functions
            auto getRowNNZ  = ampCSRParams->getRowNNZFunction();
            auto getRowCols = ampCSRParams->getRowColsFunction();

            // get NNZ counts and trigger allocations in blocks
            std::vector<lidx_t> nnz_diag( nrows ), nnz_offd( nrows );
            for ( lidx_t n = 0; n < nrows; ++n ) {
                getRowNNZ( d_first_row + n, nnz_diag[n], nnz_offd[n] );
            }
            d_diag_matrix->setNNZ( nnz_diag );
            d_offd_matrix->setNNZ( nnz_offd );


            // Fill in column indices
            for ( lidx_t n = 0; n < nrows; ++n ) {
                const auto rs_diag = d_diag_matrix->d_row_starts[n];
                auto cols_diag     = &( d_diag_matrix->d_cols[rs_diag] );
                const auto rs_offd = d_offd_matrix->d_row_starts[n];
                auto cols_offd =
                    d_offd_matrix->d_is_empty ? nullptr : &( d_offd_matrix->d_cols[rs_offd] );
                getRowCols( d_first_row + n, cols_diag, cols_offd );
            }

            // trigger re-packing of columns and convert to local cols
            globalToLocalColumns();
        }

    } else {
        AMP_ERROR( "Check supplied MatrixParameters object" );
    }

    // get total nnz count
    d_nnz = d_diag_matrix->d_nnz + d_offd_matrix->d_nnz;

    // determine if DOFManagers and CommLists need to be (re)created
    if ( d_nnz > 0 ) {
        resetDOFManagers();
    }

    d_is_square = ( d_leftDOFManager->numGlobalDOF() == d_rightDOFManager->numGlobalDOF() );
}

template<typename Policy, class Allocator, class DiagMatrixData>
CSRMatrixData<Policy, Allocator, DiagMatrixData>::~CSRMatrixData()
{
    AMPManager::decrementResource( "CSRMatrixData" );
}

template<typename Policy, class Allocator, class DiagMatrixData>
std::shared_ptr<MatrixData>
CSRMatrixData<Policy, Allocator, DiagMatrixData>::cloneMatrixData() const
{
    std::shared_ptr<CSRMatrixData> cloneData;

    cloneData = std::make_shared<CSRMatrixData<Policy, Allocator, DiagMatrixData>>();

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

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixData<Policy, Allocator, DiagMatrixData>::setNNZ( const std::vector<lidx_t> &nnz_diag,
                                                               const std::vector<lidx_t> &nnz_offd )
{
    // forward to internal blocks to get the internals allocated
    d_diag_matrix->setNNZ( nnz_diag );
    d_offd_matrix->setNNZ( nnz_offd );
    d_nnz = d_diag_matrix->d_nnz + d_offd_matrix->d_nnz;
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixData<Policy, Allocator, DiagMatrixData>::globalToLocalColumns()
{
    d_diag_matrix->globalToLocalColumns();
    d_offd_matrix->globalToLocalColumns();
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixData<Policy, Allocator, DiagMatrixData>::resetDOFManagers()
{
    auto comm = getComm();

    // There is no easy way to determine the remote DOFs and comm pattern
    // for the left vector. This side's DOFManager/CommList are rarely used
    // so we only create them if they don't exist
    bool need_left_dm = !d_leftDOFManager;
    bool need_left_cl = !d_leftCommList;
    if ( comm.anyReduce( need_left_cl ) ) {
        auto cl_params         = std::make_shared<CommunicationListParameters>();
        cl_params->d_comm      = comm;
        cl_params->d_localsize = d_last_row - d_first_row;
        d_leftCommList         = std::make_shared<CommunicationList>( cl_params );
    }
    if ( comm.anyReduce( need_left_dm ) ) {
        d_leftDOFManager =
            std::make_shared<Discretization::DOFManager>( d_last_row - d_first_row, comm );
    }

    // Right DOFManager and CommList used mre often. Replacing DOFManager has
    // poor side effects and is only done if necessary. The CommList on the other
    // hand must contain only the minimal set of remote DOFs to avoid useless
    // communication
    bool need_right_dm = !d_rightDOFManager;
    if ( d_rightDOFManager ) {
        auto dm_rdofs = d_rightDOFManager->getRemoteDOFs();
        if ( static_cast<size_t>( d_offd_matrix->numUniqueColumns() ) > dm_rdofs.size() ) {
            // too few rdofs, must remake it
            need_right_dm = true;
        } else {
            // test if needed DOFs contained in DOFManagers remotes?
        }
    }
    bool need_right_cl = !d_rightCommList;
    if ( d_rightCommList ) {
        // right CL does exist, get remote dofs and test them
        auto cl_rdofs = d_rightCommList->getGhostIDList();
        if ( static_cast<size_t>( d_offd_matrix->numUniqueColumns() ) != cl_rdofs.size() ) {
            // wrong number of rdofs, no further testing needed
            need_right_cl = true;
        } else {
            // test if individual DOFs match?
        }
    }

    if ( comm.anyReduce( need_right_cl ) ) {
        auto cl_params         = std::make_shared<CommunicationListParameters>();
        cl_params->d_comm      = comm;
        cl_params->d_localsize = d_last_col - d_first_col;
        d_offd_matrix->getColumnMap( cl_params->d_remote_DOFs );
        d_rightCommList = std::make_shared<CommunicationList>( cl_params );
    }

    if ( comm.anyReduce( need_right_dm ) ) {
        d_rightDOFManager = std::make_shared<Discretization::DOFManager>(
            d_last_col - d_first_col, comm, d_rightCommList->getGhostIDList() );
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
std::shared_ptr<DiagMatrixData> CSRMatrixData<Policy, Allocator, DiagMatrixData>::subsetRows(
    const std::vector<gidx_t> &rows ) const
{
    auto sub_matrix = std::make_shared<DiagMatrixData>( nullptr,
                                                        d_memory_location,
                                                        0,
                                                        static_cast<gidx_t>( rows.size() ),
                                                        0,
                                                        numGlobalColumns(),
                                                        true );

    // count nnz per row and write into sub matrix directly
    // also check that passed in rows are in ascending order and owned here
    gidx_t row_prev = rows[0];
    for ( size_t n = 0; n < rows.size(); ++n ) {
        AMP_DEBUG_ASSERT( n == 0 || rows[n] > row_prev );
        AMP_DEBUG_ASSERT( d_first_row <= rows[n] && rows[n] < d_last_row );
        const auto row_loc = static_cast<lidx_t>( rows[n] - d_first_row );
        sub_matrix->d_row_starts[n] +=
            ( d_diag_matrix->d_row_starts[row_loc + 1] - d_diag_matrix->d_row_starts[row_loc] );
        if ( !d_offd_matrix->d_is_empty ) {
            sub_matrix->d_row_starts[n] +=
                ( d_offd_matrix->d_row_starts[row_loc + 1] - d_offd_matrix->d_row_starts[row_loc] );
        }
        row_prev = rows[n];
    }

    // call setNNZ with accumulation on to convert counts and allocate internally
    sub_matrix->setNNZ( true );

    // Loop back over diag/offd and copy in marked rows
    for ( size_t n = 0; n < rows.size(); ++n ) {
        const auto row_loc = static_cast<lidx_t>( rows[n] - d_first_row );
        lidx_t pos         = 0;
        for ( lidx_t k = d_diag_matrix->d_row_starts[row_loc];
              k < d_diag_matrix->d_row_starts[row_loc + 1];
              ++k ) {
            sub_matrix->d_cols[k + pos] =
                d_diag_matrix->localToGlobal( d_diag_matrix->d_cols_loc[k] );
            sub_matrix->d_coeffs[k + pos] = d_diag_matrix->d_coeffs[k];
            ++pos;
        }
        if ( !d_offd_matrix->d_is_empty ) {
            for ( lidx_t k = d_offd_matrix->d_row_starts[row_loc];
                  k < d_offd_matrix->d_row_starts[row_loc + 1];
                  ++k ) {
                sub_matrix->d_cols[k + pos] =
                    d_offd_matrix->localToGlobal( d_offd_matrix->d_cols_loc[k] );
                sub_matrix->d_coeffs[k + pos] = d_offd_matrix->d_coeffs[k];
                ++pos;
            }
        }
    }

    return sub_matrix;
}

/** \brief  Extract subset of each row containing global columns in some range
 * \param[in] idx_lo  Lower global column index (inclusive)
 * \param[in] idx_up  Upper global column index (exclusive)
 * \return  shared_ptr to CSRLocalMatrixData holding the extracted nonzeros
 * \details  Returned matrix concatenates contributions for both diag and
 * offd components. Row and column extents are inherited from this matrix,
 * but are neither sorted nor converted to local indices.
 */
template<typename Policy, class Allocator, class DiagMatrixData>
std::shared_ptr<DiagMatrixData>
CSRMatrixData<Policy, Allocator, DiagMatrixData>::subsetCols( const gidx_t idx_lo,
                                                              const gidx_t idx_up ) const
{
    AMP_DEBUG_ASSERT( idx_up > idx_lo );

    auto sub_matrix = std::make_shared<DiagMatrixData>(
        nullptr, d_memory_location, d_first_row, d_last_row, d_first_col, d_last_col, true );

    // count nnz within each row that lie in the given range
    const auto nrows = static_cast<lidx_t>( d_last_row - d_first_row );
    for ( lidx_t row = 0; row < nrows; ++row ) {
        for ( lidx_t k = d_diag_matrix->d_row_starts[row]; k < d_diag_matrix->d_row_starts[row + 1];
              ++k ) {
            const auto col = d_diag_matrix->localToGlobal( d_diag_matrix->d_cols_loc[k] );
            if ( idx_lo <= col && col < idx_up ) {
                sub_matrix->d_row_starts[row]++;
            }
        }
        if ( !d_offd_matrix->d_is_empty ) {
            for ( lidx_t k = d_offd_matrix->d_row_starts[row];
                  k < d_offd_matrix->d_row_starts[row + 1];
                  ++k ) {
                const auto col = d_offd_matrix->localToGlobal( d_offd_matrix->d_cols_loc[k] );
                if ( idx_lo <= col && col < idx_up ) {
                    sub_matrix->d_row_starts[row]++;
                }
            }
        }
    }

    // call setNNZ with accumulation on to convert counts and allocate internally
    sub_matrix->setNNZ( true );

    // loop back over rows and write desired entries into sub matrix
    for ( lidx_t row = 0; row < nrows; ++row ) {
        lidx_t pos = 0;
        for ( lidx_t k = d_diag_matrix->d_row_starts[row]; k < d_diag_matrix->d_row_starts[row + 1];
              ++k ) {
            const auto col = d_diag_matrix->localToGlobal( d_diag_matrix->d_cols_loc[k] );
            if ( idx_lo <= col && col < idx_up ) {
                sub_matrix->d_cols[k + pos]   = col;
                sub_matrix->d_coeffs[k + pos] = d_diag_matrix->d_coeffs[k];
                ++pos;
            }
        }
        if ( !d_offd_matrix->d_is_empty ) {
            for ( lidx_t k = d_offd_matrix->d_row_starts[row];
                  k < d_offd_matrix->d_row_starts[row + 1];
                  ++k ) {
                const auto col = d_offd_matrix->localToGlobal( d_offd_matrix->d_cols_loc[k] );
                if ( idx_lo <= col && col < idx_up ) {
                    sub_matrix->d_cols[k + pos]   = col;
                    sub_matrix->d_coeffs[k + pos] = d_offd_matrix->d_coeffs[k];
                    ++pos;
                }
            }
        }
    }

    return sub_matrix;
}

template<typename Policy, class Allocator, class DiagMatrixData>
std::shared_ptr<MatrixData> CSRMatrixData<Policy, Allocator, DiagMatrixData>::transpose() const
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixData<Policy, Allocator, DiagMatrixData>::getRowByGlobalID(
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

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixData<Policy, Allocator, DiagMatrixData>::getValuesByGlobalID(
    size_t num_rows,
    size_t num_cols,
    size_t *rows,
    size_t *cols,
    void *vals,
    [[maybe_unused]] const typeID &id ) const
{
    AMP_DEBUG_INSIST( getTypeID<scalar_t>() == id,
                      "CSRMatrixData::getValuesByGlobalID called with inconsistent typeID" );

    AMP_DEBUG_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                      "CSRMatrixData::getValuesByGlobalID not implemented for device memory" );

    auto values = reinterpret_cast<scalar_t *>( vals );

    // zero out values
    for ( size_t i = 0; i < num_rows * num_cols; i++ ) {
        values[i] = 0.0;
    }

    // get values row-by-row from the enclosed blocks
    lidx_t start_pos = 0;
    for ( size_t nr = 0; nr < num_rows; ++nr ) {
        const auto local_row = static_cast<lidx_t>( rows[nr] - d_first_row );
        d_diag_matrix->getValuesByGlobalID(
            local_row, num_cols, &cols[start_pos], &values[start_pos] );
        d_offd_matrix->getValuesByGlobalID(
            local_row, num_cols, &cols[start_pos], &values[start_pos] );
        start_pos += num_cols;
    }
}

// The two getValues functions above can directly forward to the diag and off diag blocks
// The addValuesByGlobalID and setValuesByGlobalID functions can't do this since
// they need to also handle the other_data case
template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixData<Policy, Allocator, DiagMatrixData>::addValuesByGlobalID(
    size_t num_rows,
    size_t num_cols,
    size_t *rows,
    size_t *cols,
    void *vals,
    [[maybe_unused]] const typeID &id )
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
            d_diag_matrix->addValuesByGlobalID( num_cols, local_row, cols, lvals );
            d_offd_matrix->addValuesByGlobalID( num_cols, local_row, cols, lvals );
        } else {
            for ( size_t icol = 0; icol < num_cols; ++icol ) {
                d_other_data[rows[i]][cols[icol]] += values[num_cols * i + icol];
            }
        }
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixData<Policy, Allocator, DiagMatrixData>::setValuesByGlobalID(
    size_t num_rows,
    size_t num_cols,
    size_t *rows,
    size_t *cols,
    void *vals,
    [[maybe_unused]] const typeID &id )
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
            d_diag_matrix->setValuesByGlobalID( num_cols, local_row, cols, lvals );
            d_offd_matrix->setValuesByGlobalID( num_cols, local_row, cols, lvals );
        } else {
            for ( size_t icol = 0; icol < num_cols; ++icol ) {
                d_ghost_data[rows[i]][cols[icol]] = values[num_cols * i + icol];
            }
        }
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
std::vector<size_t>
CSRMatrixData<Policy, Allocator, DiagMatrixData>::getColumnIDs( size_t row ) const
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

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixData<Policy, Allocator, DiagMatrixData>::setOtherData(
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
                if constexpr ( std::is_same_v<gidx_t, size_t> ) {
                    addValuesByGlobalID( 1u,
                                         1u,
                                         &aggregateRows[i],
                                         &aggregateCols[i],
                                         &aggregateData[i],
                                         getTypeID<scalar_t>() );
                } else {
                    size_t row = static_cast<size_t>( aggregateRows[i] );
                    size_t col = static_cast<size_t>( aggregateCols[i] );
                    addValuesByGlobalID(
                        1u, 1u, &row, &col, &aggregateData[i], getTypeID<scalar_t>() );
                }
            }
        }
    } else {

        if ( t == AMP::LinearAlgebra::ScatterType::CONSISTENT_SET ) {
            for ( int i = 0; i != totDataLen; i++ ) {
                if ( ( aggregateRows[i] >= d_first_row ) && ( aggregateRows[i] < d_last_row ) ) {
                    if constexpr ( std::is_same_v<gidx_t, size_t> ) {
                        setValuesByGlobalID( 1u,
                                             1u,
                                             &aggregateRows[i],
                                             &aggregateCols[i],
                                             &aggregateData[i],
                                             getTypeID<scalar_t>() );
                    } else {
                        size_t row = static_cast<size_t>( aggregateRows[i] );
                        size_t col = static_cast<size_t>( aggregateCols[i] );
                        setValuesByGlobalID(
                            1u, 1u, &row, &col, &aggregateData[i], getTypeID<scalar_t>() );
                    }
                }
            }
        }
    }

    other_data.clear();
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixData<Policy, Allocator, DiagMatrixData>::makeConsistent(
    AMP::LinearAlgebra::ScatterType t )
{
    if ( t == AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD )
        setOtherData( d_other_data, AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
    else
        setOtherData( d_ghost_data, AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}

template<typename Policy, class Allocator, class DiagMatrixData>
std::shared_ptr<Discretization::DOFManager>
CSRMatrixData<Policy, Allocator, DiagMatrixData>::getRightDOFManager() const
{
    return d_rightDOFManager;
}

template<typename Policy, class Allocator, class DiagMatrixData>
std::shared_ptr<Discretization::DOFManager>
CSRMatrixData<Policy, Allocator, DiagMatrixData>::getLeftDOFManager() const
{
    return d_leftDOFManager;
}

template<typename Policy, class Allocator, class DiagMatrixData>
std::shared_ptr<CommunicationList>
CSRMatrixData<Policy, Allocator, DiagMatrixData>::getRightCommList() const
{
    return d_rightCommList;
}

template<typename Policy, class Allocator, class DiagMatrixData>
std::shared_ptr<CommunicationList>
CSRMatrixData<Policy, Allocator, DiagMatrixData>::getLeftCommList() const
{
    return d_leftCommList;
}

/********************************************************
 * Get the number of rows/columns in the matrix          *
 ********************************************************/
template<typename Policy, class Allocator, class DiagMatrixData>
size_t CSRMatrixData<Policy, Allocator, DiagMatrixData>::numLocalRows() const
{
    return static_cast<size_t>( d_last_row - d_first_row );
}

template<typename Policy, class Allocator, class DiagMatrixData>
size_t CSRMatrixData<Policy, Allocator, DiagMatrixData>::numGlobalRows() const
{
    AMP_ASSERT( d_leftDOFManager );
    return d_leftDOFManager->numGlobalDOF();
}

template<typename Policy, class Allocator, class DiagMatrixData>
size_t CSRMatrixData<Policy, Allocator, DiagMatrixData>::numLocalColumns() const
{
    return static_cast<size_t>( d_last_col - d_first_col );
}

template<typename Policy, class Allocator, class DiagMatrixData>
size_t CSRMatrixData<Policy, Allocator, DiagMatrixData>::numGlobalColumns() const
{
    AMP_ASSERT( d_rightDOFManager );
    return d_rightDOFManager->numGlobalDOF();
}

/********************************************************
 * Get iterators                                         *
 ********************************************************/
template<typename Policy, class Allocator, class DiagMatrixData>
size_t CSRMatrixData<Policy, Allocator, DiagMatrixData>::beginRow() const
{
    return static_cast<size_t>( d_first_row );
}

template<typename Policy, class Allocator, class DiagMatrixData>
size_t CSRMatrixData<Policy, Allocator, DiagMatrixData>::endRow() const
{
    return static_cast<size_t>( d_last_row );
}

template<typename Policy, class Allocator, class DiagMatrixData>
size_t CSRMatrixData<Policy, Allocator, DiagMatrixData>::beginCol() const
{
    return static_cast<size_t>( d_first_col );
}

template<typename Policy, class Allocator, class DiagMatrixData>
size_t CSRMatrixData<Policy, Allocator, DiagMatrixData>::endCol() const
{
    return static_cast<size_t>( d_last_col );
}

} // namespace AMP::LinearAlgebra

#endif
