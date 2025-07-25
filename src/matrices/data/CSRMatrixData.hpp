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

#include "ProfilerApp.h"

#include <type_traits>

namespace AMP::LinearAlgebra {

/********************************************************
 * Constructors/Destructor                              *
 ********************************************************/
template<typename Config>
CSRMatrixData<Config>::CSRMatrixData()
    : d_memory_location( AMP::Utilities::getAllocatorMemoryType<allocator_type>() )
{
    AMPManager::incrementResource( "CSRMatrixData" );
}

template<typename Config>
CSRMatrixData<Config>::CSRMatrixData( std::shared_ptr<MatrixParametersBase> params )
    : MatrixData( params ),
      d_memory_location( AMP::Utilities::getAllocatorMemoryType<allocator_type>() )
{
    PROFILE( "CSRMatrixData::CSRMatrixData" );

    AMPManager::incrementResource( "CSRMatrixData" );

    AMP_INSIST(
        d_memory_location != AMP::Utilities::MemoryType::device,
        "CSRMatrixData and CSRSerialMatrixData do not support pure-device memory locations yet" );

    // Figure out what kind of parameters object we have
    // Note: matParams always true if ampCSRParams is by inheritance
    auto rawCSRParams = std::dynamic_pointer_cast<RawCSRMatrixParameters<Config>>( d_pParameters );
    auto ampCSRParams = std::dynamic_pointer_cast<AMPCSRMatrixParameters<Config>>( d_pParameters );
    auto matParams    = std ::dynamic_pointer_cast<MatrixParameters>( d_pParameters );

    if ( rawCSRParams ) {

        // Simplest initialization, extract row/column bounds and pass through to diag/offd
        d_first_row = rawCSRParams->d_first_row;
        d_last_row  = rawCSRParams->d_last_row;
        d_first_col = rawCSRParams->d_first_col;
        d_last_col  = rawCSRParams->d_last_col;

        // Construct on/off diag blocks
        d_diag_matrix = std::make_shared<localmatrixdata_t>(
            params, d_memory_location, d_first_row, d_last_row, d_first_col, d_last_col, true );
        d_offd_matrix = std::make_shared<localmatrixdata_t>(
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
        d_diag_matrix = std::make_shared<localmatrixdata_t>(
            params, d_memory_location, d_first_row, d_last_row, d_first_col, d_last_col, true );
        d_offd_matrix = std::make_shared<localmatrixdata_t>(
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
    resetDOFManagers();

    d_is_square = ( d_leftDOFManager->numGlobalDOF() == d_rightDOFManager->numGlobalDOF() );
}

template<typename Config>
CSRMatrixData<Config>::~CSRMatrixData()
{
    AMPManager::decrementResource( "CSRMatrixData" );
}

template<typename Config>
std::shared_ptr<MatrixData> CSRMatrixData<Config>::cloneMatrixData() const
{
    std::shared_ptr<CSRMatrixData> cloneData;

    cloneData = std::make_shared<CSRMatrixData<Config>>();

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

template<typename Config>
std::shared_ptr<MatrixData> CSRMatrixData<Config>::transpose() const
{
    PROFILE( "CSRMatrixData::transpose" );

    std::shared_ptr<CSRMatrixData> transposeData;

    transposeData = std::make_shared<CSRMatrixData<Config>>();

    // copy fields from current, take care to swap L/R and rows/cols
    transposeData->d_is_square       = d_is_square;
    transposeData->d_first_row       = d_first_col;
    transposeData->d_last_row        = d_last_col;
    transposeData->d_first_col       = d_first_row;
    transposeData->d_last_col        = d_last_row;
    transposeData->d_leftDOFManager  = d_rightDOFManager;
    transposeData->d_rightDOFManager = d_leftDOFManager;
    transposeData->d_leftCommList    = d_rightCommList;
    transposeData->d_rightCommList   = d_leftCommList;

    // Parameters object is touchier, should also swap its internal L/R fields
    // Matrix can be built from many different MatrixParameters classes
    // There is no need to explicitly match the same type of parameters class
    // that was used initially
    transposeData->d_pParameters =
        std::make_shared<MatrixParameters>( d_rightDOFManager,
                                            d_leftDOFManager,
                                            getComm(),
                                            d_pParameters->getRightVariable(),
                                            d_pParameters->getLeftVariable(),
                                            d_rightCommList,
                                            d_leftCommList );

    transposeData->d_diag_matrix = d_diag_matrix->transpose( transposeData->d_pParameters );
    if ( getComm().getSize() > 1 ) {
        transposeData->d_offd_matrix = transposeOffd( transposeData->d_pParameters );
    } else {
        transposeData->d_offd_matrix =
            std::make_shared<localmatrixdata_t>( transposeData->d_pParameters,
                                                 d_memory_location,
                                                 d_first_col,
                                                 d_last_col,
                                                 d_first_row,
                                                 d_last_row,
                                                 false );
    }

    // total number of local non-zeros not the same, offd block can change
    transposeData->d_nnz =
        transposeData->d_diag_matrix->d_nnz + transposeData->d_offd_matrix->d_nnz;

    // matrix blocks will not have correct ordering within rows and still
    // have their global indices present. Call g2l to fix that.
    transposeData->globalToLocalColumns();
    transposeData->resetDOFManagers( true );

    return transposeData;
}

template<typename Config>
std::shared_ptr<CSRLocalMatrixData<Config>>
CSRMatrixData<Config>::transposeOffd( std::shared_ptr<MatrixParametersBase> params ) const
{
    PROFILE( "CSRMatrixData::transposeOffd" );

    // make a matrix communicator based on right comm list
    CSRMatrixCommunicator<Config> mat_comm( d_rightCommList );

    // extract info from offd block
    auto col_map = d_offd_matrix->getColumnMap();
    auto num_unq = d_offd_matrix->numUniqueColumns();

    // Get the partition from right comm list and test
    // which blocks need to be created
    // use the fact that partitions and col_map are sorted
    auto partition = d_rightCommList->getPartition();
    std::vector<int> dest_ranks;
    int rd = 0;
    for ( lidx_t n = 0; n < num_unq && rd < static_cast<int>( partition.size() ); ++n ) {
        const auto col  = col_map[n];
        auto part_start = static_cast<gidx_t>( rd == 0 ? 0 : partition[rd - 1] );
        auto part_end   = static_cast<gidx_t>( partition[rd] );
        if ( col < part_start ) {
            // by sorting the partition containing this index should
            // already be flagged
            continue;
        } else if ( col < part_end ) {
            // Found column in current partition, flag it and increment
            dest_ranks.push_back( rd );
            ++rd;
        } else {
            // Found column past current partition
            // increment partition until it is contained
            while ( col >= part_end ) {
                ++rd;
                AMP_DEBUG_ASSERT( rd < static_cast<int>( partition.size() ) );
                part_start = part_end;
                part_end   = partition[rd];
            }
            // insert and increment again
            dest_ranks.push_back( rd );
            ++rd;
        }
    }

    // Create blocks by subsetting on columns and send to owners
    std::map<int, std::shared_ptr<localmatrixdata_t>> send_blocks;
    for ( const auto rd : dest_ranks ) {
        const auto part_start = static_cast<gidx_t>( rd == 0 ? 0 : partition[rd - 1] );
        const auto part_end   = static_cast<gidx_t>( partition[rd] );
        auto block            = subsetCols( part_start, part_end );
        send_blocks.insert( { rd, block->transpose( params ) } );
    }
    mat_comm.sendMatrices( send_blocks );

    // receive all blocks needed here
    // swap this ranks row/col extents to get transpose's extents
    auto recv_blocks = mat_comm.recvMatrices( d_first_col, d_last_col, d_first_row, d_last_row );

    // handle edge case of no recv'd matrices (e.g. parallel matrix is block diagonal)
    if ( recv_blocks.size() == 0 ) {
        return std::make_shared<localmatrixdata_t>(
            params, d_memory_location, d_first_row, d_last_row, d_first_col, d_last_col, false );
    }

    // return horizontal concatenation of recv'd blocks
    return localmatrixdata_t::ConcatHorizontal( params, recv_blocks );
}

template<typename Config>
void CSRMatrixData<Config>::setNNZ( lidx_t tot_nnz_diag, lidx_t tot_nnz_offd )
{
    PROFILE( "CSRMatrixData::setNNZ" );

    // forward to internal blocks to get the internals allocated
    d_diag_matrix->setNNZ( tot_nnz_diag );
    d_offd_matrix->setNNZ( tot_nnz_offd );
    d_nnz = d_diag_matrix->d_nnz + d_offd_matrix->d_nnz;
}

template<typename Config>
void CSRMatrixData<Config>::setNNZ( const std::vector<lidx_t> &nnz_diag,
                                    const std::vector<lidx_t> &nnz_offd )
{
    PROFILE( "CSRMatrixData::setNNZ" );

    // forward to internal blocks to get the internals allocated
    d_diag_matrix->setNNZ( nnz_diag );
    d_offd_matrix->setNNZ( nnz_offd );
    d_nnz = d_diag_matrix->d_nnz + d_offd_matrix->d_nnz;
}

template<typename Config>
void CSRMatrixData<Config>::setNNZ( bool do_accum )
{
    PROFILE( "CSRMatrixData::setNNZ" );

    // forward to internal blocks to get the internals allocated
    d_diag_matrix->setNNZ( do_accum );
    d_offd_matrix->setNNZ( do_accum );
    d_nnz = d_diag_matrix->d_nnz + d_offd_matrix->d_nnz;
}

template<typename Config>
void CSRMatrixData<Config>::globalToLocalColumns()
{
    PROFILE( "CSRMatrixData::globalToLocalColumns" );

    d_diag_matrix->globalToLocalColumns();
    d_offd_matrix->globalToLocalColumns();
}

template<typename Config>
void CSRMatrixData<Config>::resetDOFManagers( bool force_right )
{
    PROFILE( "CSRMatrixData::resetDOFManagers" );

    auto comm = getComm();

    // There is no easy way to determine the remote DOFs and comm pattern
    // for the left vector. This side's DOFManager/CommList are rarely used
    // so we only create them if they don't exist
    if ( !d_leftCommList ) {
        auto cl_params         = std::make_shared<CommunicationListParameters>();
        cl_params->d_comm      = comm;
        cl_params->d_localsize = d_last_row - d_first_row;
        d_leftCommList         = std::make_shared<CommunicationList>( cl_params );
    }
    if ( !d_leftDOFManager ) {
        d_leftDOFManager =
            std::make_shared<Discretization::DOFManager>( d_last_row - d_first_row, comm );
    }

    // Right DOFManager and CommList used mre often. Replacing DOFManager has
    // poor side effects and is only done if necessary. The CommList on the other
    // hand must contain only the minimal set of remote DOFs to avoid useless
    // communication
    bool need_right_dm = !d_rightDOFManager || force_right;
    if ( d_rightDOFManager ) {
        auto dm_rdofs = d_rightDOFManager->getRemoteDOFs();
        if ( static_cast<size_t>( d_offd_matrix->numUniqueColumns() ) > dm_rdofs.size() ) {
            // too few rdofs, must remake it
            need_right_dm = true;
        } else {
            // test if needed DOFs contained in DOFManagers remotes?
        }
    }
    bool need_right_cl = !d_rightCommList || force_right;
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

template<typename Config>
std::shared_ptr<CSRLocalMatrixData<Config>>
CSRMatrixData<Config>::subsetRows( const std::vector<gidx_t> &rows ) const
{
    PROFILE( "CSRMatrixData::subsetRows" );

    auto sub_matrix = std::make_shared<localmatrixdata_t>( nullptr,
                                                           d_memory_location,
                                                           0,
                                                           static_cast<gidx_t>( rows.size() ),
                                                           0,
                                                           numGlobalColumns(),
                                                           true );

    // count nnz per row and write into sub matrix directly
    // also check that passed in rows are in ascending order and owned here
    [[maybe_unused]] gidx_t row_prev = rows[0];
    for ( size_t n = 0; n < rows.size(); ++n ) {
        AMP_DEBUG_ASSERT( n == 0 || rows[n] > row_prev );
        AMP_DEBUG_ASSERT( d_first_row <= rows[n] && rows[n] < d_last_row );
        const auto row_loc = static_cast<lidx_t>( rows[n] - d_first_row );
        sub_matrix->d_row_starts[n] =
            ( d_diag_matrix->d_row_starts[row_loc + 1] - d_diag_matrix->d_row_starts[row_loc] );
        if ( !d_offd_matrix->d_is_empty ) {
            sub_matrix->d_row_starts[n] +=
                ( d_offd_matrix->d_row_starts[row_loc + 1] - d_offd_matrix->d_row_starts[row_loc] );
        }
        row_prev = rows[n];
    }

    // call setNNZ with accumulation on to convert counts and allocate internally
    sub_matrix->setNNZ( true );
    if ( sub_matrix->d_nnz == 0 ) {
        std::cout << "Got zero nnz" << std::endl;
        return sub_matrix;
    }

    // Loop back over diag/offd and copy in marked rows
    for ( size_t n = 0; n < rows.size(); ++n ) {
        const auto row_loc = static_cast<lidx_t>( rows[n] - d_first_row );
        lidx_t pos         = sub_matrix->d_row_starts[n];
        for ( lidx_t k = d_diag_matrix->d_row_starts[row_loc];
              k < d_diag_matrix->d_row_starts[row_loc + 1];
              ++k ) {
            sub_matrix->d_cols[pos] = d_diag_matrix->localToGlobal( d_diag_matrix->d_cols_loc[k] );
            sub_matrix->d_coeffs[pos] = d_diag_matrix->d_coeffs[k];
            ++pos;
        }
        if ( !d_offd_matrix->d_is_empty ) {
            for ( lidx_t k = d_offd_matrix->d_row_starts[row_loc];
                  k < d_offd_matrix->d_row_starts[row_loc + 1];
                  ++k ) {
                sub_matrix->d_cols[pos] =
                    d_offd_matrix->localToGlobal( d_offd_matrix->d_cols_loc[k] );
                sub_matrix->d_coeffs[pos] = d_offd_matrix->d_coeffs[k];
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
template<typename Config>
std::shared_ptr<CSRLocalMatrixData<Config>>
CSRMatrixData<Config>::subsetCols( const gidx_t idx_lo, const gidx_t idx_up ) const
{
    PROFILE( "CSRMatrixData::subsetCols" );

    AMP_DEBUG_ASSERT( idx_up > idx_lo );

    auto sub_matrix = std::make_shared<localmatrixdata_t>(
        nullptr, d_memory_location, d_first_row, d_last_row, idx_lo, idx_up, true );

    // count nnz within each row that lie in the given range
    const auto nrows = static_cast<lidx_t>( d_last_row - d_first_row );
    for ( lidx_t row = 0; row < nrows; ++row ) {
        sub_matrix->d_row_starts[row] = 0;
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
        lidx_t pos = sub_matrix->d_row_starts[row];
        for ( lidx_t k = d_diag_matrix->d_row_starts[row]; k < d_diag_matrix->d_row_starts[row + 1];
              ++k ) {
            const auto col = d_diag_matrix->localToGlobal( d_diag_matrix->d_cols_loc[k] );
            if ( idx_lo <= col && col < idx_up ) {
                sub_matrix->d_cols[pos]   = col;
                sub_matrix->d_coeffs[pos] = d_diag_matrix->d_coeffs[k];
                ++pos;
            }
        }
        if ( !d_offd_matrix->d_is_empty ) {
            for ( lidx_t k = d_offd_matrix->d_row_starts[row];
                  k < d_offd_matrix->d_row_starts[row + 1];
                  ++k ) {
                const auto col = d_offd_matrix->localToGlobal( d_offd_matrix->d_cols_loc[k] );
                if ( idx_lo <= col && col < idx_up ) {
                    sub_matrix->d_cols[pos]   = col;
                    sub_matrix->d_coeffs[pos] = d_offd_matrix->d_coeffs[k];
                    ++pos;
                }
            }
        }
        AMP_DEBUG_ASSERT( pos == sub_matrix->d_row_starts[row + 1] );
    }

    return sub_matrix;
}

template<typename Config>
void CSRMatrixData<Config>::getRowByGlobalID( size_t row,
                                              std::vector<size_t> &cols,
                                              std::vector<double> &vals ) const
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

template<typename Config>
void CSRMatrixData<Config>::getValuesByGlobalID( size_t num_rows,
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
template<typename Config>
void CSRMatrixData<Config>::addValuesByGlobalID( size_t num_rows,
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

template<typename Config>
void CSRMatrixData<Config>::setValuesByGlobalID( size_t num_rows,
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

template<typename Config>
std::vector<size_t> CSRMatrixData<Config>::getColumnIDs( size_t row ) const
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

template<typename Config>
void CSRMatrixData<Config>::setOtherData( std::map<gidx_t, std::map<gidx_t, scalar_t>> &other_data,
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

template<typename Config>
void CSRMatrixData<Config>::makeConsistent( AMP::LinearAlgebra::ScatterType t )
{
    if ( t == AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD )
        setOtherData( d_other_data, AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
    else
        setOtherData( d_ghost_data, AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}

template<typename Config>
std::shared_ptr<Discretization::DOFManager> CSRMatrixData<Config>::getRightDOFManager() const
{
    return d_rightDOFManager;
}

template<typename Config>
std::shared_ptr<Discretization::DOFManager> CSRMatrixData<Config>::getLeftDOFManager() const
{
    return d_leftDOFManager;
}

template<typename Config>
std::shared_ptr<CommunicationList> CSRMatrixData<Config>::getRightCommList() const
{
    return d_rightCommList;
}

template<typename Config>
std::shared_ptr<CommunicationList> CSRMatrixData<Config>::getLeftCommList() const
{
    return d_leftCommList;
}

/********************************************************
 * Get the number of rows/columns in the matrix          *
 ********************************************************/
template<typename Config>
size_t CSRMatrixData<Config>::numLocalRows() const
{
    return static_cast<size_t>( d_last_row - d_first_row );
}

template<typename Config>
size_t CSRMatrixData<Config>::numGlobalRows() const
{
    AMP_ASSERT( d_leftDOFManager );
    return d_leftDOFManager->numGlobalDOF();
}

template<typename Config>
size_t CSRMatrixData<Config>::numLocalColumns() const
{
    return static_cast<size_t>( d_last_col - d_first_col );
}

template<typename Config>
size_t CSRMatrixData<Config>::numGlobalColumns() const
{
    AMP_ASSERT( d_rightDOFManager );
    return d_rightDOFManager->numGlobalDOF();
}

/********************************************************
 * Get iterators                                         *
 ********************************************************/
template<typename Config>
size_t CSRMatrixData<Config>::beginRow() const
{
    return static_cast<size_t>( d_first_row );
}

template<typename Config>
size_t CSRMatrixData<Config>::endRow() const
{
    return static_cast<size_t>( d_last_row );
}

template<typename Config>
size_t CSRMatrixData<Config>::beginCol() const
{
    return static_cast<size_t>( d_first_col );
}

template<typename Config>
size_t CSRMatrixData<Config>::endCol() const
{
    return static_cast<size_t>( d_last_col );
}

} // namespace AMP::LinearAlgebra

#endif
