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

#if defined( USE_CUDA ) || defined( USE_HIP )
    #include <thrust/device_vector.h>
    #include <thrust/execution_policy.h>
    #include <thrust/for_each.h>
    #include <thrust/inner_product.h>
#endif


#include <algorithm>
#include <iterator>
#include <memory>
#include <numeric>
#include <set>
#include <type_traits>

namespace AMP::LinearAlgebra {

/********************************************************
 * Constructor/Destructor helper functions              *
 ********************************************************/
template<typename Policy>
bool isColValid( typename Policy::gidx_t col,
                 bool is_diag,
                 typename Policy::gidx_t first_col,
                 typename Policy::gidx_t last_col )
{
    bool dValid  = is_diag && ( first_col <= col && col < last_col );
    bool odValid = !is_diag && ( col < first_col || last_col <= col );
    return ( dValid || odValid );
}

/********************************************************
 * Constructors/Destructor                              *
 ********************************************************/
template<typename Policy, class Allocator>
CSRMatrixData<Policy, Allocator>::CSRMatrixData()
{
    AMPManager::incrementResource( "CSRMatrixData" );
}

template<typename Policy, class Allocator>
CSRMatrixData<Policy, Allocator>::CSRSerialMatrixData::CSRSerialMatrixData(
    const CSRMatrixData<Policy, Allocator> &outer )
    : d_outer( outer )
{
    AMPManager::incrementResource( "CSRSerialMatrixData" );
}

template<typename Policy, class Allocator>
CSRMatrixData<Policy, Allocator>::CSRMatrixData( std::shared_ptr<MatrixParametersBase> params )
    : MatrixData( params )
{

    AMPManager::incrementResource( "CSRMatrixData" );
    auto csrParams = std::dynamic_pointer_cast<CSRMatrixParameters<Policy>>( d_pParameters );
    auto matParams = std ::dynamic_pointer_cast<MatrixParameters>( d_pParameters );

    d_memory_location = d_pParameters->d_memory_location;

    if ( csrParams ) {

        // add check for memory location etc and migrate if necessary
        d_is_square = csrParams->d_is_square;
        d_first_row = csrParams->d_first_row;
        d_last_row  = csrParams->d_last_row;
        d_first_col = csrParams->d_first_col;
        d_last_col  = csrParams->d_last_col;

        size_t N = d_last_row - d_first_row;

        if ( d_memory_location != AMP::Utilities::MemoryType::device ) {
            // Construct on/off diag blocks
            d_diag_matrix     = std::make_shared<CSRSerialMatrixData>( *this, params, true );
            d_off_diag_matrix = std::make_shared<CSRSerialMatrixData>( *this, params, false );

            // get total nnz count
            d_nnz = d_diag_matrix->d_nnz + d_off_diag_matrix->d_nnz;

            // collect off-diagonal entries and create right dof manager
            std::vector<size_t> remote_dofs;
            for ( lidx_t i = 0; i < d_off_diag_matrix->d_nnz; ++i ) {
                remote_dofs.push_back( d_off_diag_matrix->d_cols[i] );
            }
            AMP::Utilities::unique( remote_dofs );
            const auto &comm = getComm();
            d_rightDOFManager =
                std::make_shared<AMP::Discretization::DOFManager>( N, comm, remote_dofs );

            if ( d_is_square ) {
                d_leftDOFManager = d_rightDOFManager;
            } else {
                AMP_ERROR( "Non-square matrices not handled at present" );
            }
        } else {
            AMP_WARNING( "CSRMatrixData: device memory handling has not been implemented yet" );
        }
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

        // send params forward to the on/off diagonal blocks
        d_diag_matrix     = std::make_shared<CSRSerialMatrixData>( *this, params, true );
        d_off_diag_matrix = std::make_shared<CSRSerialMatrixData>( *this, params, false );
        d_nnz             = d_diag_matrix->d_nnz + d_off_diag_matrix->d_nnz;

    } else {
        AMP_ERROR( "Check supplied MatrixParameter object" );
    }
}

template<typename Policy, class Allocator>
CSRMatrixData<Policy, Allocator>::~CSRMatrixData()
{
    AMPManager::decrementResource( "CSRMatrixData" );
}

/********************************************************
 * Constructors/Destructor for nested class             *
 ********************************************************/
template<typename Policy, class Allocator>
CSRMatrixData<Policy, Allocator>::CSRSerialMatrixData::CSRSerialMatrixData(
    const CSRMatrixData<Policy, Allocator> &outer,
    std::shared_ptr<MatrixParametersBase> params,
    bool is_diag )
    : d_outer( outer )
{
    AMPManager::incrementResource( "CSRSerialMatrixData" );
    d_pParameters  = params;
    auto csrParams = std::dynamic_pointer_cast<CSRMatrixParameters<Policy>>( d_pParameters );
    auto matParams = std ::dynamic_pointer_cast<MatrixParameters>( d_pParameters );

    d_memory_location = d_pParameters->d_memory_location;
    d_is_diag         = is_diag;

    // Number of rows owned by this rank
    d_num_rows = outer.d_last_row - outer.d_first_row;

    if ( csrParams ) {
        // Pull out block specific parameters
        auto &blParams = d_is_diag ? csrParams->d_diag : csrParams->d_off_diag;

        // memory not managed here regardless of block type (except row starts)
        d_own_data = false;

        // copy in data pointers
        d_nnz_per_row = blParams.d_nnz_per_row;
        d_row_starts  = blParams.d_row_starts;
        d_cols        = blParams.d_cols;
        d_cols_loc    = blParams.d_cols_loc;
        d_coeffs      = blParams.d_coeffs;
        d_nnz_pad     = d_is_diag ? 0 : csrParams->d_nnz_pad;

        // count nnz and decide if block is empty
        d_nnz      = std::accumulate( d_nnz_per_row, d_nnz_per_row + d_num_rows, 0 );
        d_is_empty = ( d_nnz == 0 );

    } else if ( matParams ) {

        // for now all matrix parameter data is assumed to be on host

        auto leftDOFManager  = matParams->getLeftDOFManager();
        auto rightDOFManager = matParams->getRightDOFManager();
        AMP_ASSERT( leftDOFManager && rightDOFManager );
        AMP_ASSERT( matParams->d_CommListLeft && matParams->d_CommListRight );

        d_is_empty = false;

        const auto &getRow = matParams->getRowFunction();
        AMP_INSIST( getRow,
                    "Explicitly defined getRow function must be present in MatrixParameters"
                    " to construct CSRMatrixData and CSRSerialMatrixData" );

        // Count number of nonzeros depending on block type
        // also track un-referenced columns if off-diagonal
        std::vector<gidx_t> colPad;
        std::set<gidx_t> colSet;
        d_nnz_pad = 0;
        d_nnz     = 0;
        for ( gidx_t i = outer.d_first_row; i < outer.d_last_row; ++i ) {
            for ( auto &&col : getRow( i ) ) {
                if ( isColValid<Policy>( col, d_is_diag, outer.d_first_col, outer.d_last_col ) ) {
                    ++d_nnz;
                    if ( !d_is_diag ) {
                        colSet.insert( col );
                    }
                }
            }
        }

        // attempt to insert all remote dofs into colSet to see which are un-referenced
        if ( !d_is_diag ) {
            auto remoteDOFs = rightDOFManager->getRemoteDOFs();
            for ( auto &&rdof : remoteDOFs ) {
                auto cs = colSet.insert( rdof );
                if ( cs.second ) {
                    // insertion success means this DOF is un-referenced
                    // add it to the padding list
                    colPad.push_back( rdof );
                    ++d_nnz;
                    ++d_nnz_pad;
                }
            }
        }

        // bail out for degenerate case with no nnz
        // may happen in off-diagonal blocks
        if ( d_nnz == 0 ) {
            d_is_empty = true;
            d_own_data = false;
            return;
        }

        // allocate internal arrays
        d_own_data    = true;
        d_nnz_per_row = lidxAllocator.allocate( d_num_rows );
        d_row_starts  = lidxAllocator.allocate( d_num_rows + 1 );
        d_cols        = gidxAllocator.allocate( d_nnz );
        d_cols_loc    = lidxAllocator.allocate( d_nnz );
        d_coeffs      = scalarAllocator.allocate( d_nnz );

        // Fill cols and nnz based on local row extents and on/off diag status
        lidx_t cli       = 0; // index into local array of columns as it is filled in
        lidx_t nnzFilled = 0;
        for ( lidx_t i = 0; i < d_num_rows; ++i ) {
            d_nnz_per_row[i] = 0;
            auto cols        = getRow( outer.d_first_row + i );
            for ( auto &&col : cols ) {
                if ( isColValid<Policy>( col, d_is_diag, outer.d_first_col, outer.d_last_col ) ) {
                    d_nnz_per_row[i]++;
                    d_cols[cli] = col;
                    if ( d_is_diag ) {
                        d_cols_loc[cli] = static_cast<lidx_t>( col - outer.d_first_col );
                    } else {
                        d_cols_loc[cli] = static_cast<lidx_t>(
                            matParams->d_CommListRight->getLocalGhostID( col ) );
                    }
                    d_coeffs[cli] = 0.0;
                    ++cli;
                    ++nnzFilled;
                }
            }
        }

        // If off-diag pad in the un-referenced ghosts to the final row
        if ( !d_is_diag ) {
            d_nnz_per_row[d_num_rows - 1] += d_nnz_pad;

            for ( auto col : colPad ) {
                d_cols[cli] = col;
                d_cols_loc[cli] =
                    static_cast<lidx_t>( matParams->d_CommListRight->getLocalGhostID( col ) );
                d_coeffs[cli] = 0.0;
                ++cli;
                ++nnzFilled;
            }
        }

        // scan nnz counts to get starting index of each row
        std::exclusive_scan( d_nnz_per_row, d_nnz_per_row + d_num_rows, d_row_starts, 0 );
        d_row_starts[d_num_rows] = d_row_starts[d_num_rows - 1] + d_nnz_per_row[d_num_rows - 1];

        // Ensure that the right number of nnz were actually filled in
        AMP_DEBUG_ASSERT( nnzFilled == d_nnz );
    } else {
        AMP_ERROR( "Check supplied MatrixParameter object" );
    }
}
template<typename Policy, class Allocator>
void CSRMatrixData<Policy, Allocator>::CSRSerialMatrixData::findColumnMap()
{
    if ( d_ncols_unq > 0 ) {
        // return if it already known
        return;
    }

    // Otherwise allocate and fill the map
    // Number of unique (global) columns is largest value in local cols
    d_ncols_unq = *( std::max_element( d_cols_loc, d_cols_loc + d_nnz ) );
    ++d_ncols_unq; // plus one for zero-based indexing

    // Map is not allocated by default
    d_cols_unq = gidxAllocator.allocate( d_ncols_unq );

    // Fill by writing in d_cols indexed by d_cols_loc
    for ( lidx_t n = 0; n < d_nnz; ++n ) {
        d_cols_unq[d_cols_loc[n]] = d_cols[n];
    }
}

template<typename Policy, class Allocator>
CSRMatrixData<Policy, Allocator>::CSRSerialMatrixData::~CSRSerialMatrixData()
{
    AMPManager::decrementResource( "CSRSerialMatrixData" );

    // Always attempt deletion of column map since it is created
    // lazily and not loaned to other instances
    gidxAllocator.deallocate( d_cols_unq, d_ncols_unq );

    // Deallocate remaining data only if this instance owns it
    if ( d_own_data ) {
        lidxAllocator.deallocate( d_row_starts, d_num_rows + 1 );
        gidxAllocator.deallocate( d_cols, d_nnz );
        lidxAllocator.deallocate( d_cols_loc, d_nnz );
        lidxAllocator.deallocate( d_nnz_per_row, d_num_rows );
        scalarAllocator.deallocate( d_coeffs, d_nnz );
    }
}

template<typename Policy, class Allocator>
std::shared_ptr<MatrixData> CSRMatrixData<Policy, Allocator>::cloneMatrixData() const
{
    std::shared_ptr<CSRMatrixData> cloneData;

    cloneData = std::make_shared<CSRMatrixData<Policy, Allocator>>();

    cloneData->d_memory_location = d_memory_location;
    cloneData->d_is_square       = d_is_square;
    cloneData->d_first_row       = d_first_row;
    cloneData->d_last_row        = d_last_row;
    cloneData->d_first_col       = d_first_col;
    cloneData->d_last_col        = d_last_col;
    cloneData->d_nnz             = d_nnz;
    cloneData->d_leftDOFManager  = d_leftDOFManager;
    cloneData->d_rightDOFManager = d_rightDOFManager;
    cloneData->d_pParameters     = d_pParameters;

    cloneData->d_diag_matrix     = d_diag_matrix->cloneMatrixData( *cloneData );
    cloneData->d_off_diag_matrix = d_off_diag_matrix->cloneMatrixData( *cloneData );

    return cloneData;
}

template<typename Policy, class Allocator>
std::shared_ptr<typename CSRMatrixData<Policy, Allocator>::CSRSerialMatrixData>
CSRMatrixData<Policy, Allocator>::CSRSerialMatrixData::cloneMatrixData(
    const CSRMatrixData<Policy, Allocator> &outer )
{
    std::shared_ptr<CSRSerialMatrixData> cloneData;

    cloneData = std::make_shared<CSRSerialMatrixData>( outer );

    cloneData->d_is_diag         = d_is_diag;
    cloneData->d_is_empty        = d_is_empty;
    cloneData->d_num_rows        = d_num_rows;
    cloneData->d_nnz             = d_nnz;
    cloneData->d_memory_location = d_memory_location;
    cloneData->d_pParameters     = d_pParameters;

    if ( !d_is_empty ) {
        cloneData->d_own_data    = true;
        cloneData->d_nnz_per_row = lidxAllocator.allocate( d_num_rows );
        cloneData->d_row_starts  = lidxAllocator.allocate( d_num_rows + 1 );
        cloneData->d_cols        = gidxAllocator.allocate( d_nnz );
        cloneData->d_cols_loc    = lidxAllocator.allocate( d_nnz );
        cloneData->d_coeffs      = scalarAllocator.allocate( d_nnz );

        if ( d_memory_location < AMP::Utilities::MemoryType::device ) {
            std::copy( d_nnz_per_row, d_nnz_per_row + d_num_rows, cloneData->d_nnz_per_row );
            std::copy( d_row_starts, d_row_starts + d_num_rows + 1, cloneData->d_row_starts );
            std::copy( d_cols, d_cols + d_nnz, cloneData->d_cols );
            std::copy( d_cols_loc, d_cols_loc + d_nnz, cloneData->d_cols_loc );
            // need to zero out coeffs so that padded region has valid data
#warning May remove fill here when padding is removed
            std::fill( d_coeffs, d_coeffs + d_nnz, 0.0 );
        } else {
#if defined( USE_CUDA ) || defined( USE_HIP )
            // I hope this is temporary. Generally, I advocate for CSRMatrixData being
            // execution space agnostic. I have some ideas for that. (Brian Romero)
            thrust::copy_n( thrust::device, d_nnz_per_row, d_num_rows, cloneData->d_nnz_per_row );
            thrust::copy_n( thrust::device, d_num_rows + d_row_starts, 1, cloneData->d_row_starts );
            thrust::copy_n( thrust::device, d_cols, d_nnz, cloneData->d_cols );
            thrust::copy_n( thrust::device, d_cols_loc, d_nnz, cloneData->d_cols_loc );
                // need to zero out coeffs so that padded region has valid data
    #warning May remove fill here when padding is removed
            thrust::fill_n( thrust::device, d_coeffs, d_nnz, 0.0 );
#else
            AMP_ERROR( "No device found!" );
#endif
        }
    } else {
        cloneData->d_own_data    = false;
        cloneData->d_nnz_per_row = nullptr;
        cloneData->d_row_starts  = nullptr;
        cloneData->d_cols        = nullptr;
        cloneData->d_cols_loc    = nullptr;
        cloneData->d_coeffs      = nullptr;
    }

    return cloneData;
}

template<typename Policy, class Allocator>
std::shared_ptr<MatrixData> CSRMatrixData<Policy, Allocator>::transpose() const
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy, class Allocator>
void CSRMatrixData<Policy, Allocator>::extractDiagonal( std::shared_ptr<Vector> buf ) const
{
    AMP_ASSERT( buf && buf->numberOfDataBlocks() == 1 ); // temporary constraint
    AMP_ASSERT( buf->isType<scalar_t>( 0 ) );

    auto *rawVecData = buf->getRawDataBlock<scalar_t>();
    auto memType     = AMP::Utilities::getMemoryType( rawVecData );
    if ( memType < AMP::Utilities::MemoryType::device ) {

        const size_t N = d_last_row - d_first_row;
        for ( size_t i = 0; i < N; ++i ) {
            const auto start = d_diag_matrix->d_row_starts[i];
            const auto end   = d_diag_matrix->d_row_starts[i + 1];
            // colums are unordered at present
            for ( lidx_t j = start; j < end; ++j ) {
                if ( d_diag_matrix->d_cols[j] == static_cast<gidx_t>( d_first_col + i ) ) {
                    rawVecData[i] = d_diag_matrix->d_coeffs[j];
                    break;
                }
            }
        }
    } else {
        AMP_ERROR(
            "CSRSerialMatrixData<Policy>::extractDiagonal not implemented for vec and matrix in "
            "different memory spaces" );
    }
}
template<typename Policy, class Allocator>
void CSRMatrixData<Policy, Allocator>::getRowByGlobalID( size_t row,
                                                         std::vector<size_t> &cols,
                                                         std::vector<double> &vals ) const
{
    AMP_INSIST( row >= static_cast<size_t>( d_first_row ) &&
                    row < static_cast<size_t>( d_last_row ),
                "row must be owned by rank" );

    auto local_row = row - d_first_row;

    // Get portion of row from diagonal matrix
    d_diag_matrix->getRowByGlobalID( local_row, cols, vals );

    // Get portion from off diagonal and append
    std::vector<size_t> od_cols;
    std::vector<double> od_vals;
    d_off_diag_matrix->getRowByGlobalID( local_row, od_cols, od_vals );
    cols.insert( cols.end(), od_cols.begin(), od_cols.end() );
    vals.insert( vals.end(), od_vals.begin(), od_vals.end() );
}

template<typename Policy, class Allocator>
void CSRMatrixData<Policy, Allocator>::getValuesByGlobalID( size_t num_rows,
                                                            size_t num_cols,
                                                            size_t *rows,
                                                            size_t *cols,
                                                            void *values,
                                                            const typeID &id ) const
{
    if ( getTypeID<scalar_t>() == id ) {
        if ( d_memory_location < AMP::Utilities::MemoryType::device ) {

            if ( num_rows == 1 && num_cols == 1 ) {

                const auto local_row = rows[0] - d_first_row;
                // Forward to internal matrices, nothing will happen if not found
                d_diag_matrix->getValuesByGlobalID( local_row, cols[0], values, id );
                d_off_diag_matrix->getValuesByGlobalID( local_row, cols[0], values, id );
            } else {
                AMP_ERROR(
                    "CSRSerialMatrixData::getValuesByGlobalID not implemented for num_rows>1 || "
                    "num_cols > 1" );
            }

        } else {
            AMP_ERROR(
                "CSRSerialMatrixData::getValuesByGlobalID not implemented for device memory" );
        }
    } else {
        AMP_ERROR( "Not implemented" );
    }
}

// The two getValues functions above can directly forward to the diag and off diag blocks
// The addValuesByGlobalID and setValuesByGlobalID functions can't do this since
// they need to also handle the other_data case
template<typename Policy, class Allocator>
void CSRMatrixData<Policy, Allocator>::addValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, void *vals, const typeID &id )
{
    if ( getTypeID<scalar_t>() != id ) {
        AMP_ERROR( "Conversion not implemented" );
    }

    if ( d_memory_location < AMP::Utilities::MemoryType::device ) {

        auto values = reinterpret_cast<const scalar_t *>( vals );

        for ( size_t i = 0u; i != num_rows; i++ ) {
            if ( rows[i] >= static_cast<size_t>( d_first_row ) &&
                 rows[i] < static_cast<size_t>( d_last_row ) ) {

                // Forward single row to diag and off diag blocks
                // auto lcols = &cols[num_cols * i];
                const auto local_row = rows[i] - d_first_row;
                auto lvals           = &values[num_cols * i];
                d_diag_matrix->addValuesByGlobalID( num_cols, local_row, cols, lvals, id );
                d_off_diag_matrix->addValuesByGlobalID( num_cols, local_row, cols, lvals, id );
            } else {
                for ( size_t icol = 0; icol < num_cols; ++icol ) {
                    d_other_data[rows[i]][cols[icol]] += values[num_cols * i + icol];
                }
            }
        }

    } else {
        AMP_ERROR( "CSRMatrixData::addValuesByGlobalID not implemented for device memory" );
    }
}

template<typename Policy, class Allocator>
void CSRMatrixData<Policy, Allocator>::setValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, void *vals, const typeID &id )
{
    if ( getTypeID<scalar_t>() != id ) {
        AMP_ERROR( "Conversion not implemented" );
    }

    if ( d_memory_location < AMP::Utilities::MemoryType::device ) {

        auto values = reinterpret_cast<const scalar_t *>( vals );

        for ( size_t i = 0u; i != num_rows; i++ ) {

            if ( rows[i] >= static_cast<size_t>( d_first_row ) &&
                 rows[i] < static_cast<size_t>( d_last_row ) ) {

                // Forward single row to diag and off diag blocks
                // auto lcols = &cols[num_cols * i];
                const auto local_row = rows[i] - d_first_row;
                auto lvals           = &values[num_cols * i];
                d_diag_matrix->setValuesByGlobalID( num_cols, local_row, cols, lvals, id );
                d_off_diag_matrix->setValuesByGlobalID( num_cols, local_row, cols, lvals, id );
            } else {
                for ( size_t icol = 0; icol < num_cols; ++icol ) {
                    d_ghost_data[rows[i]][cols[icol]] = values[num_cols * i + icol];
                }
            }
        }

    } else {
        AMP_ERROR( "CSRMatrixData::addValuesByGlobalID not implemented for device memory" );
    }
}

template<typename Policy, class Allocator>
std::vector<size_t> CSRMatrixData<Policy, Allocator>::getColumnIDs( size_t row ) const
{
    AMP_INSIST( row >= static_cast<size_t>( d_first_row ) &&
                    row < static_cast<size_t>( d_last_row ),
                "row must be owned by rank" );
    AMP_INSIST( d_diag_matrix, "diag matrix must exist" );
    auto local_row              = row - d_first_row;
    std::vector<size_t> cols    = d_diag_matrix->getColumnIDs( local_row );
    std::vector<size_t> od_cols = d_off_diag_matrix->getColumnIDs( local_row );
    cols.insert( cols.end(), od_cols.begin(), od_cols.end() );
    return cols;
}

template<typename Policy, class Allocator>
void CSRMatrixData<Policy, Allocator>::CSRSerialMatrixData::getRowByGlobalID(
    const size_t local_row, std::vector<size_t> &cols, std::vector<double> &values ) const
{
    // Don't do anything on empty matrices
    if ( d_is_empty ) {
        return;
    }

    if ( d_memory_location < AMP::Utilities::MemoryType::device ) {
        const size_t last_row = d_num_rows - 1;
        const auto row_offset = static_cast<size_t>( local_row );
        const auto offset     = std::accumulate( d_nnz_per_row, d_nnz_per_row + row_offset, 0 );
        auto n                = d_nnz_per_row[row_offset];
        if ( local_row == last_row ) {
            n -= d_nnz_pad;
        }

        cols.resize( n );
        values.resize( n );

        if constexpr ( std::is_same_v<size_t, gidx_t> ) {
            std::copy( &d_cols[offset], &d_cols[offset] + n, cols.begin() );
        } else {
            std::transform( &d_cols[offset],
                            &d_cols[offset] + n,
                            cols.begin(),
                            []( size_t col ) -> gidx_t { return col; } );
        }

        if constexpr ( std::is_same_v<double, scalar_t> ) {
            std::copy( &d_coeffs[offset], &d_coeffs[offset] + n, values.begin() );
        } else {
            std::transform( &d_coeffs[offset],
                            &d_coeffs[offset] + n,
                            values.begin(),
                            []( size_t val ) -> scalar_t { return val; } );
        }
    } else {
        AMP_ERROR( "CSRSerialMatrixData::getRowByGlobalID not implemented for device memory" );
    }
}

template<typename Policy, class Allocator>
void CSRMatrixData<Policy, Allocator>::CSRSerialMatrixData::getValuesByGlobalID(
    const size_t local_row, const size_t col, void *values, const typeID &id ) const
{
    // Don't do anything on empty matrices
    if ( d_is_empty ) {
        return;
    }

    if ( getTypeID<scalar_t>() != id ) {
        AMP_ERROR( "Conversion not implemented" );
    }

    const size_t last_row = d_num_rows - 1;
    const auto start      = d_row_starts[local_row];
    auto end              = d_row_starts[local_row + 1];
    if ( local_row == last_row ) {
#warning This code is jank. Go fix the constructor that needs padding.
        end -= d_nnz_pad;
    }

    for ( lidx_t i = start; i < end; ++i ) {
        if ( d_cols[i] == static_cast<gidx_t>( col ) ) {
            *( reinterpret_cast<scalar_t *>( values ) ) = d_coeffs[i];
        }
    }
}

template<typename Policy, class Allocator>
void CSRMatrixData<Policy, Allocator>::CSRSerialMatrixData::addValuesByGlobalID(
    const size_t num_cols,
    const size_t local_row,
    const size_t *cols,
    const scalar_t *vals,
    const typeID &id )
{
    if ( d_is_empty ) {
        return;
    }

    if ( getTypeID<scalar_t>() != id ) {
        AMP_ERROR( "Conversion not implemented" );
    }

    const size_t last_row = d_num_rows - 1;
    const auto start      = d_row_starts[local_row];
    auto end              = d_row_starts[local_row + 1];
    if ( local_row == last_row ) {
#warning This code is jank. Go fix the constructor that needs padding.
        end -= d_nnz_pad;
    }

    // Inefficient because we don't assume order
    // not sure it's worth optimizing for our use cases
    for ( size_t icol = 0; icol < num_cols; ++icol ) {
        for ( lidx_t j = start; j < end; ++j ) {
            if ( d_cols[j] == static_cast<gidx_t>( cols[icol] ) ) {
                d_coeffs[j] += vals[icol];
            }
        }
    }
}

template<typename Policy, class Allocator>
void CSRMatrixData<Policy, Allocator>::CSRSerialMatrixData::setValuesByGlobalID(
    const size_t num_cols,
    const size_t local_row,
    const size_t *cols,
    const scalar_t *vals,
    const typeID &id )
{
    if ( d_is_empty ) {
        return;
    }

    if ( getTypeID<scalar_t>() != id ) {
        AMP_ERROR( "Conversion not implemented" );
    }

    const size_t last_row = d_num_rows - 1;
    const auto start      = d_row_starts[local_row];
    auto end              = d_row_starts[local_row + 1];
    if ( local_row == last_row ) {
#warning This code is jank. Go fix the constructor that needs padding.
        end -= d_nnz_pad;
    }

    // Inefficient because we don't assume order
    // not sure it's worth optimizing for our use cases
    for ( size_t icol = 0; icol < num_cols; ++icol ) {
        for ( lidx_t j = start; j < end; ++j ) {
            if ( d_cols[j] == static_cast<gidx_t>( cols[icol] ) ) {
                d_coeffs[j] = vals[icol];
                if ( j > ( d_nnz - d_nnz_pad ) ) {
                    AMP_INSIST( d_coeffs[j] == 0.0, " Assigning non-zero to padded location" );
                }
            }
        }
    }
}

template<typename Policy, class Allocator>
std::vector<size_t>
CSRMatrixData<Policy, Allocator>::CSRSerialMatrixData::getColumnIDs( const size_t local_row ) const
{
    // Don't do anything on empty matrices
    if ( d_is_empty ) {
        return std::vector<size_t>();
    }

    AMP_INSIST( d_cols && d_nnz_per_row, "Must be initialized" );

    if ( d_memory_location < AMP::Utilities::MemoryType::device ) {

        std::vector<size_t> cols;
        const size_t last_row = d_num_rows - 1;
        const auto row_offset = static_cast<size_t>( local_row );
        const auto offset     = d_row_starts[local_row];
        auto n                = d_nnz_per_row[row_offset];

        if ( local_row == last_row ) {
#warning This code is jank. Go fix the constructor that needs padding.
            n -= d_nnz_pad;
        }

        if constexpr ( std::is_same_v<size_t, gidx_t> ) {
            std::copy( &d_cols[offset], &d_cols[offset] + n, std::back_inserter( cols ) );
        } else {
            std::transform( &d_cols[offset],
                            &d_cols[offset] + n,
                            std::back_inserter( cols ),
                            []( size_t col ) -> gidx_t { return col; } );
        }
        return cols;
    } else {
        AMP_ERROR( "CSRSerialMatrixData:getColumnIDs not implemented for device memory" );
    }
}

template<typename Policy, class Allocator>
void CSRMatrixData<Policy, Allocator>::setOtherData(
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

template<typename Policy, class Allocator>
void CSRMatrixData<Policy, Allocator>::makeConsistent( AMP::LinearAlgebra::ScatterType t )
{
    if ( t == AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD )
        setOtherData( d_other_data, AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
    else
        setOtherData( d_ghost_data, AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}

template<typename Policy, class Allocator>
std::shared_ptr<Discretization::DOFManager>
CSRMatrixData<Policy, Allocator>::getRightDOFManager() const
{
    return d_rightDOFManager;
}

template<typename Policy, class Allocator>
std::shared_ptr<Discretization::DOFManager>
CSRMatrixData<Policy, Allocator>::getLeftDOFManager() const
{
    return d_leftDOFManager;
}

/********************************************************
 * Get the number of rows/columns in the matrix          *
 ********************************************************/
template<typename Policy, class Allocator>
size_t CSRMatrixData<Policy, Allocator>::numLocalRows() const
{
    return static_cast<size_t>( d_last_row - d_first_row );
}

template<typename Policy, class Allocator>
size_t CSRMatrixData<Policy, Allocator>::numGlobalRows() const
{
    AMP_ASSERT( d_leftDOFManager );
    return d_leftDOFManager->numGlobalDOF();
}

template<typename Policy, class Allocator>
size_t CSRMatrixData<Policy, Allocator>::numLocalColumns() const
{
    return static_cast<size_t>( d_last_col - d_first_col );
}

template<typename Policy, class Allocator>
size_t CSRMatrixData<Policy, Allocator>::numGlobalColumns() const
{
    AMP_ASSERT( d_rightDOFManager );
    return d_rightDOFManager->numGlobalDOF();
}

/********************************************************
 * Get iterators                                         *
 ********************************************************/
template<typename Policy, class Allocator>
size_t CSRMatrixData<Policy, Allocator>::beginRow() const
{
    return static_cast<size_t>( d_first_row );
}

template<typename Policy, class Allocator>
size_t CSRMatrixData<Policy, Allocator>::endRow() const
{
    return static_cast<size_t>( d_last_row );
}

} // namespace AMP::LinearAlgebra

#endif
