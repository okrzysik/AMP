#ifndef included_AMP_CSRLocalMatrixData_hpp
#define included_AMP_CSRLocalMatrixData_hpp

#include "AMP/AMP_TPLs.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/AMPCSRMatrixParameters.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/matrices/RawCSRMatrixParameters.h"
#include "AMP/matrices/data/CSRLocalMatrixData.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Algorithms.h"
#include "AMP/utils/Utilities.h"

#include <numeric>
#include <set>
#include <type_traits>

namespace AMP::LinearAlgebra {

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

template<typename size_type, class data_allocator>
std::shared_ptr<typename data_allocator::value_type[]> sharedArrayBuilder( size_type N,
                                                                           data_allocator &alloc )
{
    AMP_DEBUG_ASSERT( std::is_integral_v<size_type> );
    return std::shared_ptr<typename data_allocator::value_type[]>(
        alloc.allocate( N ), [N, &alloc]( auto p ) -> void { alloc.deallocate( p, N ); } );
}

template<typename data_type>
std::shared_ptr<data_type[]> sharedArrayWrapper( data_type *raw_array )
{
    return std::shared_ptr<data_type[]>( raw_array, []( auto p ) -> void { (void) p; } );
}

template<typename Policy, class Allocator>
CSRLocalMatrixData<Policy, Allocator>::CSRLocalMatrixData(
    std::shared_ptr<MatrixParametersBase> params,
    AMP::Utilities::MemoryType memory_location,
    typename Policy::gidx_t first_row,
    typename Policy::gidx_t last_row,
    typename Policy::gidx_t first_col,
    typename Policy::gidx_t last_col,
    bool is_diag )
    : d_memory_location( memory_location ),
      d_first_row( first_row ),
      d_last_row( last_row ),
      d_first_col( first_col ),
      d_last_col( last_col ),
      d_is_diag( is_diag ),
      d_num_rows( last_row - first_row )
{
    AMPManager::incrementResource( "CSRLocalMatrixData" );

    // Figure out what kind of parameters object we have
    // Note: matParams always true if ampCSRParams is by inheritance
    auto rawCSRParams = std::dynamic_pointer_cast<RawCSRMatrixParameters<Policy>>( params );
    auto ampCSRParams = std::dynamic_pointer_cast<AMPCSRMatrixParameters<Policy>>( params );
    auto matParams    = std ::dynamic_pointer_cast<MatrixParameters>( params );

    if ( rawCSRParams ) {
        // Pull out block specific parameters
        auto &blParams = d_is_diag ? rawCSRParams->d_diag : rawCSRParams->d_off_diag;

        if ( blParams.d_row_starts == nullptr ) {
            d_is_empty = true;
            return;
        }

        // count nnz and decide if block is empty
        d_nnz = blParams.d_row_starts[d_num_rows];

        if ( d_nnz == 0 ) {
            d_is_empty = true;
            return;
        }
        d_is_empty = false;

        // Wrap raw pointers from blParams to match internal
        // shared_ptr<T[]> type
        d_row_starts = sharedArrayWrapper( blParams.d_row_starts );
        d_cols       = sharedArrayWrapper( blParams.d_cols );
        d_coeffs     = sharedArrayWrapper( blParams.d_coeffs );
    } else if ( matParams ) {
        // Getting device memory support in this constructor mode will be very challenging
        AMP_ASSERT( d_memory_location != AMP::Utilities::MemoryType::device );

        // can always allocate row starts without external information
        d_row_starts = sharedArrayBuilder( d_num_rows + 1, d_lidxAllocator );
        AMP::Utilities::Algorithms<lidx_t>::fill_n( d_row_starts.get(), d_num_rows + 1, 0 );

        const auto &getRow = matParams->getRowFunction();

        if ( !getRow || ampCSRParams ) {
            // Initialization not desired or not possible
            // can be set later by calling setNNZ and filling d_cols in some fashion
            d_nnz      = 0;
            d_is_empty = true;
            return;
        }

        // Count number of nonzeros per row and total
        d_nnz = 0;
        for ( gidx_t i = d_first_row; i < d_last_row; ++i ) {
            lidx_t valid_nnz = 0;
            for ( auto &&col : getRow( i ) ) {
                if ( isColValid<Policy>( col, d_is_diag, d_first_col, d_last_col ) ) {
                    ++valid_nnz;
                }
            }
            d_row_starts[i - d_first_row] = valid_nnz;
            d_nnz += valid_nnz;
        }

        // bail out for degenerate case with no nnz
        // may happen in off-diagonal blocks
        if ( d_nnz == 0 ) {
            d_is_empty = true;
            return;
        }
        d_is_empty = false;

        // Allocate internal arrays
        d_cols   = sharedArrayBuilder( d_nnz, d_gidxAllocator );
        d_coeffs = sharedArrayBuilder( d_nnz, d_scalarAllocator );

        // Fill cols and nnz based on local row extents and on/off diag status
        lidx_t nnzFilled = 0;
        lidx_t nnzCached = d_row_starts[0];
        d_row_starts[0]  = 0;
        for ( lidx_t row = 0; row < d_num_rows; ++row ) {
            // do exclusive scan on nnz per row along the way
            lidx_t rs             = d_row_starts[row] + nnzCached;
            nnzCached             = d_row_starts[row + 1];
            d_row_starts[row + 1] = rs;
            // fill in valid columns from getRow function
            auto cols = getRow( d_first_row + row );
            for ( auto &&col : cols ) {
                if ( isColValid<Policy>( col, d_is_diag, d_first_col, d_last_col ) ) {
                    d_cols[nnzFilled]   = col;
                    d_coeffs[nnzFilled] = 0.0;
                    ++nnzFilled;
                }
            }
        }

        // Ensure that the right number of nnz were actually filled in
        AMP_DEBUG_ASSERT( nnzFilled == d_nnz );
    } else {
        // The parameters object is allowed to be null
        // In this case matrices will stay purely local (e.g. not parts of
        // an enclosing CSRMatrixData object). This is used for remote blocks
        // in SpGEMM
        d_nnz        = 0;
        d_is_empty   = true;
        d_row_starts = sharedArrayBuilder( d_num_rows + 1, d_lidxAllocator );
        AMP::Utilities::Algorithms<lidx_t>::fill_n( d_row_starts.get(), d_num_rows + 1, 0 );
        return;
    }

    // fill in local column indices
    globalToLocalColumns();
}

template<typename Policy, class Allocator>
std::shared_ptr<CSRLocalMatrixData<Policy, Allocator>>
CSRLocalMatrixData<Policy, Allocator>::ConcatHorizontal(
    std::shared_ptr<MatrixParametersBase> params,
    std::map<int, std::shared_ptr<CSRLocalMatrixData<Policy, Allocator>>> blocks )
{
    PROFILE( "CSRLocalMatrixData::ConcatHorizontal" );

    AMP_INSIST( blocks.size() > 0, "Attempted to concatenate empty set of blocks" );

    // Verify that all have matching row/col starts/stops
    // Blocks must have valid global columns present
    // Count total number of non-zeros in each row from combination.
    auto block           = ( *blocks.begin() ).second;
    const auto mem_loc   = block->d_memory_location;
    const auto first_row = block->d_first_row;
    const auto last_row  = block->d_last_row;
    const auto nrows     = static_cast<lidx_t>( last_row - first_row );
    const auto first_col = block->d_first_col;
    const auto last_col  = block->d_last_col;
    std::vector<lidx_t> row_nnz( last_row - first_row, 0 );
    for ( auto it : blocks ) {
        block = it.second;
        AMP_INSIST( first_row == block->d_first_row && last_row == block->d_last_row &&
                        first_col == block->d_first_col && last_col == block->d_last_col,
                    "Blocks to concatenate must have compatible layouts" );
        AMP_INSIST( block->d_cols.get(), "Blocks to concatenate must have global columns" );
        AMP_INSIST( mem_loc == block->d_memory_location,
                    "Blocks to concatenate must be in same memory space" );
        for ( lidx_t row = 0; row < nrows; ++row ) {
            row_nnz[row] += ( block->d_row_starts[row + 1] - block->d_row_starts[row] );
        }
    }

    // Create empty matrix and trigger allocations to match
    auto concat_matrix = std::make_shared<CSRLocalMatrixData<Policy, Allocator>>(
        params, mem_loc, first_row, last_row, first_col, last_col, false );
    concat_matrix->setNNZ( row_nnz );

    // set row_nnz back to zeros to use as counters while appending entries
    std::fill( row_nnz.begin(), row_nnz.end(), 0 );

    // loop back over blocks and write into new matrix
    for ( auto it : blocks ) {
        block = it.second;
        for ( lidx_t row = 0; row < nrows; ++row ) {
            for ( auto n = block->d_row_starts[row]; n < block->d_row_starts[row + 1]; ++n ) {
                const auto rs                     = concat_matrix->d_row_starts[row];
                const auto ctr                    = row_nnz[row];
                concat_matrix->d_cols[rs + ctr]   = block->d_cols[n];
                concat_matrix->d_coeffs[rs + ctr] = block->d_coeffs[n];
                row_nnz[row]++;
            }
        }
    }

    return concat_matrix;
}

template<typename Policy, class Allocator>
std::shared_ptr<CSRLocalMatrixData<Policy, Allocator>>
CSRLocalMatrixData<Policy, Allocator>::ConcatVertical(
    std::shared_ptr<MatrixParametersBase> params,
    std::map<int, std::shared_ptr<CSRLocalMatrixData<Policy, Allocator>>> blocks,
    const gidx_t first_col,
    const gidx_t last_col,
    const bool is_diag )
{
    PROFILE( "CSRLocalMatrixData::ConcatVertical" );

    AMP_INSIST( blocks.size() > 0, "Attempted to concatenate empty set of blocks" );

    // count number of rows and check compatibility of blocks
    auto block         = ( *blocks.begin() ).second;
    const auto mem_loc = block->d_memory_location;
    lidx_t num_rows    = 0;
    bool all_empty     = block->isEmpty();
    for ( auto it : blocks ) {
        block = it.second;
        AMP_DEBUG_INSIST( mem_loc == block->d_memory_location,
                          "Blocks to concatenate must be in same memory space" );
        num_rows += block->d_num_rows;
        all_empty = all_empty && block->isEmpty();
    }

    // extreme edge case where every requested row happened to be empty
    if ( all_empty ) {
        return nullptr;
    }

    // create output matrix
    auto concat_matrix = std::make_shared<CSRLocalMatrixData<Policy, Allocator>>(
        params, mem_loc, 0, num_rows, first_col, last_col, is_diag );

    // Count total number of non-zeros in each row from combination.
    lidx_t cat_row = 0; // counter for which row we are on in concat_matrix
    for ( auto it : blocks ) {
        block = it.second;
        // loop over rows and count NZs that fall in/out of column range
        for ( lidx_t brow = 0; brow < block->d_num_rows; ++brow ) {
            lidx_t row_nnz = 0;
            for ( lidx_t k = block->d_row_starts[brow]; k < block->d_row_starts[brow + 1]; ++k ) {
                const auto c      = block->d_cols[k];
                const bool inside = first_col <= c && c < last_col;
                const bool valid  = ( is_diag && inside ) || ( !is_diag && !inside );
                if ( valid ) {
                    ++row_nnz;
                }
            }
            concat_matrix->d_row_starts[cat_row] = row_nnz;
            ++cat_row;
        }
    }

    // Trigger allocations
    concat_matrix->setNNZ( true );

    // loop over blocks again and write into new matrix
    cat_row = 0;
    for ( auto it : blocks ) {
        block = it.second;
        if ( !block->d_is_empty ) {
            for ( lidx_t brow = 0; brow < block->d_num_rows; ++brow ) {
                lidx_t cat_pos = concat_matrix->d_row_starts[cat_row];
                for ( auto k = block->d_row_starts[brow]; k < block->d_row_starts[brow + 1]; ++k ) {
                    const auto c      = block->d_cols[k];
                    const bool inside = first_col <= c && c < last_col;
                    const bool valid  = ( is_diag && inside ) || ( !is_diag && !inside );
                    if ( valid ) {
                        concat_matrix->d_cols[cat_pos]   = c;
                        concat_matrix->d_coeffs[cat_pos] = block->d_coeffs[k];
                        ++cat_pos;
                    }
                }
                ++cat_row;
            }
        }
    }

    return concat_matrix;
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::swapDataFields(
    CSRLocalMatrixData<Policy, Allocator> &other )
{
    // swap metadata
    const auto o_is_empty  = other.d_is_empty;
    const auto o_nnz       = other.d_nnz;
    const auto o_ncols_unq = other.d_ncols_unq;
    other.d_is_empty       = d_is_empty;
    other.d_nnz            = d_nnz;
    other.d_ncols_unq      = d_ncols_unq;
    d_is_empty             = o_is_empty;
    d_nnz                  = o_nnz;
    d_ncols_unq            = o_ncols_unq;
    // swap fields
    d_row_starts.swap( other.d_row_starts );
    d_cols.swap( other.d_cols );
    d_cols_loc.swap( other.d_cols_loc );
    d_cols_unq.swap( other.d_cols_unq );
    d_coeffs.swap( other.d_coeffs );
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::globalToLocalColumns()
{
    PROFILE( "CSRLocalMatrixData::globalToLocalColumns" );

    if ( d_is_empty || d_cols.get() == nullptr ) {
        // gToL either trivially not needed or has already been called
        return;
    }

    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRLocalMatrixData::globalToLocalColumns not implemented on device yet" );

    // Columns easier to sort before converting to local
    // and defining unq cols in offd easier if globals are sorted
    sortColumns();

    // local columns always owned internally
    d_cols_loc = sharedArrayBuilder( d_nnz, d_lidxAllocator );

    if ( d_is_diag ) {
        for ( lidx_t n = 0; n < d_nnz; ++n ) {
            d_cols_loc[n] = static_cast<lidx_t>( d_cols[n] - d_first_col );
        }
    } else {
        // for offd setup column map as part of the process

        // first make a copy of the global columns and sort them
        // as a whole. This is different from the sortColumns call
        // above that acts within a row, where this jumbles rows
        // together.
        auto cols_tmp = sharedArrayBuilder( d_nnz, d_gidxAllocator );
        AMP::Utilities::Algorithms<gidx_t>::copy_n( d_cols.get(), d_nnz, cols_tmp.get() );
        std::sort( cols_tmp.get(), cols_tmp.get() + d_nnz );

        // count unique entries, allocate uniques, fill uniques
        auto col_curr = cols_tmp[0];
        lidx_t pos    = 1;
        for ( lidx_t n = 1; n < d_nnz; ++n ) {
            if ( cols_tmp[n] != col_curr ) {
                col_curr = cols_tmp[n];
                ++pos;
            }
        }
        d_ncols_unq   = pos;
        d_cols_unq    = sharedArrayBuilder( d_ncols_unq, d_gidxAllocator );
        col_curr      = cols_tmp[0];
        d_cols_unq[0] = col_curr;
        pos           = 1;
        for ( lidx_t n = 1; n < d_nnz; ++n ) {
            if ( cols_tmp[n] != col_curr ) {
                col_curr        = cols_tmp[n];
                d_cols_unq[pos] = col_curr;
                ++pos;
            }
        }
        cols_tmp.reset();

        // copy and modify from AMP::Utilities::findfirst to suit task
        const gidx_t *cols_unq = d_cols_unq.get();
        const lidx_t ncols_unq = d_ncols_unq;
        auto bsearch           = [cols_unq, ncols_unq]( gidx_t gc ) -> lidx_t {
            AMP_DEBUG_ASSERT( cols_unq[0] <= gc && gc <= cols_unq[ncols_unq - 1] );
            lidx_t lower = 0, upper = ncols_unq - 1, idx;
            while ( ( upper - lower ) > 1 ) {
                idx = ( upper + lower ) / 2;
                if ( cols_unq[idx] == gc ) {
                    return idx;
                } else if ( cols_unq[idx] > gc ) {
                    upper = idx - 1;
                } else {
                    lower = idx + 1;
                }
            }
            return gc == cols_unq[upper] ? upper : lower;
        };

        // find all local column indices from set
        for ( lidx_t n = 0; n < d_nnz; ++n ) {
            d_cols_loc[n] = bsearch( d_cols[n] );
            AMP_DEBUG_ASSERT( d_cols_loc[n] < d_ncols_unq );
            AMP_DEBUG_ASSERT( d_cols_unq[d_cols_loc[n]] == d_cols[n] );
        }
    }

    // free global cols as they should not be used from here on out
    d_cols.reset();
}

template<typename Policy, class Allocator>
typename Policy::gidx_t
CSRLocalMatrixData<Policy, Allocator>::localToGlobal( const typename Policy::lidx_t loc_id ) const
{
    if ( d_is_diag ) {
        return static_cast<typename Policy::gidx_t>( loc_id ) + d_first_col;
    } else {
        return d_cols_unq[loc_id];
    }
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::sortColumns()
{
    PROFILE( "CSRLocalMatrixData::sortColumns" );

    typedef std::tuple<gidx_t, scalar_t> tuple_t;

    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRSerialMatrixData::sortColumns not implemented for device memory" );

    if ( d_is_empty ) {
        return;
    }

    AMP_DEBUG_INSIST( d_cols.get() != nullptr,
                      "CSRLocalMatrixData::sortColumns Access to global columns required" );

    std::vector<tuple_t> rTpl;
    for ( lidx_t row = 0; row < d_num_rows; ++row ) {
        const auto rs      = d_row_starts[row];
        const auto row_len = d_row_starts[row + 1] - rs;
        if ( row_len == 0 )
            continue;

        // enlarge temp vector of tuples if needed
        if ( row_len > static_cast<lidx_t>( rTpl.size() ) ) {
            rTpl.resize( row_len );
        }

        // pack local column and coeff into array of tuples
        for ( lidx_t k = 0; k < row_len; ++k ) {
            rTpl[k] = std::make_tuple( d_cols[rs + k], d_coeffs[rs + k] );
        }

        // slightly different sorting criteria for on and off diagonal blocks
        if ( d_is_diag ) {
            const gidx_t diag_idx = d_first_col + static_cast<gidx_t>( row );
            // diag block puts diag entry first, then ascending order on local col
            std::sort( rTpl.data(),
                       rTpl.data() + row_len,
                       [diag_idx]( const tuple_t &a, const tuple_t &b ) -> bool {
                           const gidx_t gca = std::get<0>( a ), gcb = std::get<0>( b );
                           return diag_idx != gcb && ( gca < gcb || gca == diag_idx );
                       } );
        } else {
            // offd block is plain ascending order on local col
            std::sort( rTpl.data(),
                       rTpl.data() + row_len,
                       []( const tuple_t &a, const tuple_t &b ) -> bool {
                           return std::get<0>( a ) < std::get<0>( b );
                       } );
        }

        // unpack now sorted array of tuples
        for ( lidx_t k = 0; k < row_len; ++k ) {
            d_cols[rs + k]   = std::get<0>( rTpl[k] );
            d_coeffs[rs + k] = std::get<1>( rTpl[k] );
        }
    }
}

template<typename Policy, class Allocator>
CSRLocalMatrixData<Policy, Allocator>::~CSRLocalMatrixData()
{
    AMPManager::decrementResource( "CSRLocalMatrixData" );
}

template<typename Policy, class Allocator>
std::shared_ptr<CSRLocalMatrixData<Policy, Allocator>>
CSRLocalMatrixData<Policy, Allocator>::cloneMatrixData()
{
    std::shared_ptr<CSRLocalMatrixData> cloneData;

    cloneData = std::make_shared<CSRLocalMatrixData>(
        nullptr, d_memory_location, d_first_row, d_last_row, d_first_col, d_last_col, d_is_diag );

    cloneData->d_is_empty = d_is_empty;
    cloneData->d_nnz      = d_nnz;

    cloneData->d_cols     = nullptr;
    cloneData->d_cols_loc = nullptr;
    cloneData->d_coeffs   = nullptr;

    if ( !d_is_empty ) {
        cloneData->d_cols_loc = sharedArrayBuilder( d_nnz, d_lidxAllocator );
        cloneData->d_coeffs   = sharedArrayBuilder( d_nnz, d_scalarAllocator );

        AMP::Utilities::Algorithms<lidx_t>::copy_n(
            d_row_starts.get(), d_num_rows + 1, cloneData->d_row_starts.get() );
        AMP::Utilities::Algorithms<lidx_t>::copy_n(
            d_cols_loc.get(), d_nnz, cloneData->d_cols_loc.get() );
        AMP::Utilities::Algorithms<scalar_t>::fill_n( cloneData->d_coeffs.get(), d_nnz, 0.0 );

        if ( d_cols.get() != nullptr ) {
            cloneData->d_cols = sharedArrayBuilder( d_nnz, d_gidxAllocator );
            AMP::Utilities::Algorithms<gidx_t>::copy_n(
                d_cols.get(), d_nnz, cloneData->d_cols.get() );
        }
        if ( d_cols_unq.get() != nullptr ) {
            cloneData->d_ncols_unq = d_ncols_unq;
            cloneData->d_cols_unq  = sharedArrayBuilder( d_ncols_unq, d_gidxAllocator );
            AMP::Utilities::Algorithms<gidx_t>::copy_n(
                d_cols_unq.get(), d_ncols_unq, cloneData->d_cols_unq.get() );
        }
    }

    return cloneData;
}

template<typename Policy, class Allocator>
std::shared_ptr<CSRLocalMatrixData<Policy, Allocator>>
CSRLocalMatrixData<Policy, Allocator>::transpose(
    std::shared_ptr<MatrixParametersBase> params ) const
{
    // create new data, note swapped rows and cols
    auto transposeData = std::make_shared<CSRLocalMatrixData>(
        params, d_memory_location, d_first_col, d_last_col, d_first_row, d_last_row, d_is_diag );

    // handle rare edge case of empty diagonal block
    if ( d_is_empty ) {
        return transposeData;
    }

    auto trans_row = [is_diag   = d_is_diag,
                      first_col = d_first_col,
                      cols      = d_cols,
                      cols_loc  = d_cols_loc,
                      cols_unq  = d_cols_unq]( const lidx_t c ) -> lidx_t {
        gidx_t col_g = 0;
        if ( cols.get() ) {
            col_g = cols[c];
        } else if ( is_diag ) {
            return cols_loc[c];
        } else {
            col_g = cols_unq[cols_loc[c]];
        }
        return col_g - first_col;
    };

    // count nnz per column and store in transpose's rowstarts array
    for ( lidx_t row = 0; row < d_num_rows; ++row ) {
        for ( lidx_t c = d_row_starts[row]; c < d_row_starts[row + 1]; ++c ) {
            const auto t_row = trans_row( c );
            transposeData->d_row_starts[t_row]++;
        }
    }

    transposeData->setNNZ( true );

    // count nnz per column again and append into each row of transpose
    // create temporary vector of counters to hold position in each row
    std::vector<lidx_t> row_ctr( transposeData->d_num_rows, 0 );
    for ( lidx_t row = 0; row < d_num_rows; ++row ) {
        for ( lidx_t c = d_row_starts[row]; c < d_row_starts[row + 1]; ++c ) {
            const auto t_row = trans_row( c );
            const auto pos   = transposeData->d_row_starts[t_row] + row_ctr[t_row];
            // local transpose only fills global cols and coeffs
            // caller responsible for creation of local columns if desired
            transposeData->d_cols[pos]   = static_cast<gidx_t>( row ) + d_first_row;
            transposeData->d_coeffs[pos] = d_coeffs[c];
            row_ctr[t_row]++;
        }
    }

    return transposeData;
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::setNNZ( lidx_t tot_nnz )
{
    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRLocalMatrixData::setNNZ not implemented on device yet" );

    d_nnz = tot_nnz;

    if ( d_nnz == 0 ) {
        d_is_empty = true;
        // nothing to do, block stays empty
        return;
    }

    // allocate and fill remaining arrays
    d_is_empty = false;
    d_cols     = sharedArrayBuilder( d_nnz, d_gidxAllocator );
    d_cols_loc = sharedArrayBuilder( d_nnz, d_lidxAllocator );
    d_coeffs   = sharedArrayBuilder( d_nnz, d_scalarAllocator );

    std::fill( d_cols.get(), d_cols.get() + d_nnz, 0 );
    std::fill( d_cols_loc.get(), d_cols_loc.get() + d_nnz, 0 );
    std::fill( d_coeffs.get(), d_coeffs.get() + d_nnz, static_cast<scalar_t>( 0.0 ) );
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::setNNZ( bool do_accum )
{
    if ( do_accum ) {
        std::exclusive_scan(
            d_row_starts.get(), d_row_starts.get() + d_num_rows + 1, d_row_starts.get(), 0 );
    }

    // total nnz in all rows of block is last entry
    d_nnz = d_row_starts[d_num_rows];

    if ( d_nnz == 0 ) {
        d_is_empty = true;
        // nothing to do, block stays empty
        return;
    }

    // allocate and fill remaining arrays
    d_is_empty = false;
    d_cols     = sharedArrayBuilder( d_nnz, d_gidxAllocator );
    d_cols_loc = sharedArrayBuilder( d_nnz, d_lidxAllocator );
    d_coeffs   = sharedArrayBuilder( d_nnz, d_scalarAllocator );

    std::fill( d_cols.get(), d_cols.get() + d_nnz, 0 );
    std::fill( d_cols_loc.get(), d_cols_loc.get() + d_nnz, 0 );
    std::fill( d_coeffs.get(), d_coeffs.get() + d_nnz, static_cast<scalar_t>( 0.0 ) );
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::setNNZ( const std::vector<lidx_t> &nnz )
{
    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRLocalMatrixData::setNNZ not implemented on device yet" );

    // copy passed nnz vector into row_starts and call internal setNNZ
    std::copy( nnz.begin(), nnz.end(), d_row_starts.get() );
    setNNZ( true );
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::getColPtrs( std::vector<gidx_t *> &col_ptrs )
{
    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRLocalMatrixData::setNNZ not implemented on device yet" );

    if ( !d_is_empty ) {
        for ( lidx_t row = 0; row < d_num_rows; ++row ) {
            col_ptrs[row] = &d_cols[d_row_starts[row]];
        }
    } else {
        for ( lidx_t row = 0; row < d_num_rows; ++row ) {
            col_ptrs[row] = nullptr;
        }
    }
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::getRowByGlobalID( const size_t local_row,
                                                              std::vector<size_t> &cols,
                                                              std::vector<double> &values ) const
{
    // Don't do anything on empty matrices
    if ( d_is_empty ) {
        return;
    }

    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRLocalMatrixData::getRowByGlobalID not implemented for device memory" );

    const auto start = d_row_starts[local_row];
    auto end         = d_row_starts[local_row + 1];
    cols.resize( end - start );
    values.resize( end - start );

    // don't store global ids, need to generate on the fly
    if ( d_is_diag ) {
        const auto first_col = d_first_col;
        std::transform( &d_cols_loc[start],
                        &d_cols_loc[end],
                        cols.begin(),
                        [&]( lidx_t col ) -> size_t { return col + first_col; } );
    } else {
        const auto cols_unq = d_cols_unq;
        std::transform( &d_cols_loc[start],
                        &d_cols_loc[end],
                        cols.begin(),
                        [&]( lidx_t col ) -> size_t { return cols_unq[col]; } );
    }

    if constexpr ( std::is_same_v<double, scalar_t> ) {
        std::copy( &d_coeffs[start], &d_coeffs[end], values.begin() );
    } else {
        std::transform( &d_coeffs[start],
                        &d_coeffs[end],
                        values.begin(),
                        []( scalar_t val ) -> double { return val; } );
    }
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::getValuesByGlobalID( const size_t local_row,
                                                                 const size_t num_cols,
                                                                 size_t *cols,
                                                                 scalar_t *values ) const
{
    // Don't do anything on empty matrices
    if ( d_is_empty ) {
        return;
    }

    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRLocalMatrixData::getValuesByGlobalID not implemented for device memory" );

    const auto start = d_row_starts[local_row];
    auto end         = d_row_starts[local_row + 1];


    for ( size_t nc = 0; nc < num_cols; ++nc ) {
        auto query_col = cols[nc];
        for ( lidx_t i = start; i < end; ++i ) {
            auto icol = d_is_diag ? ( d_first_col + d_cols_loc[i] ) : ( d_cols_unq[d_cols_loc[i]] );
            if ( icol == static_cast<gidx_t>( query_col ) ) {
                values[nc] = d_coeffs[i];
            }
        }
    }
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::addValuesByGlobalID( const size_t num_cols,
                                                                 const size_t local_row,
                                                                 const size_t *cols,
                                                                 const scalar_t *vals )
{
    if ( d_is_empty ) {
        return;
    }

    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRLocalMatrixData::addValuesByGlobalID not implemented for device memory" );

    const auto start = d_row_starts[local_row];
    auto end         = d_row_starts[local_row + 1];

    if ( d_is_diag ) {
        for ( size_t icol = 0; icol < num_cols; ++icol ) {
            for ( lidx_t j = start; j < end; ++j ) {
                if ( d_first_col + d_cols_loc[j] == static_cast<gidx_t>( cols[icol] ) ) {
                    d_coeffs[j] += vals[icol];
                    break;
                }
            }
        }
    } else {
        for ( size_t icol = 0; icol < num_cols; ++icol ) {
            for ( lidx_t j = start; j < end; ++j ) {
                if ( d_cols_unq[d_cols_loc[j]] == static_cast<gidx_t>( cols[icol] ) ) {
                    d_coeffs[j] += vals[icol];
                    break;
                }
            }
        }
    }
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::setValuesByGlobalID( const size_t num_cols,
                                                                 const size_t local_row,
                                                                 const size_t *cols,
                                                                 const scalar_t *vals )
{
    if ( d_is_empty ) {
        return;
    }

    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRLocalMatrixData::setValuesByGlobalID not implemented for device memory" );

    const auto start = d_row_starts[local_row];
    auto end         = d_row_starts[local_row + 1];

    if ( d_is_diag ) {
        for ( size_t icol = 0; icol < num_cols; ++icol ) {
            for ( lidx_t j = start; j < end; ++j ) {
                if ( d_first_col + d_cols_loc[j] == static_cast<gidx_t>( cols[icol] ) ) {
                    d_coeffs[j] = vals[icol];
                    break;
                }
            }
        }
    } else {
        for ( size_t icol = 0; icol < num_cols; ++icol ) {
            for ( lidx_t j = start; j < end; ++j ) {
                if ( d_cols_unq[d_cols_loc[j]] == static_cast<gidx_t>( cols[icol] ) ) {
                    d_coeffs[j] = vals[icol];
                    break;
                }
            }
        }
    }
}

template<typename Policy, class Allocator>
std::vector<size_t>
CSRLocalMatrixData<Policy, Allocator>::getColumnIDs( const size_t local_row ) const
{
    // Don't do anything on empty matrices
    if ( d_is_empty ) {
        return std::vector<size_t>();
    }

    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRLocalMatrixData::getColumnIDs not implemented for device memory" );

    AMP_INSIST( d_cols_loc && d_row_starts,
                "CSRLocalMatrixData::getColumnIDs nnz layout must be initialized" );

    const auto start = d_row_starts[local_row];
    const auto end   = d_row_starts[local_row + 1];
    std::vector<size_t> cols( end - start, 0 );

    // don't store global ids, need to generate on the fly
    if ( d_is_diag ) {
        const auto first_col = d_first_col;
        std::transform( &d_cols_loc[start],
                        &d_cols_loc[end],
                        cols.begin(),
                        [&]( lidx_t col ) -> size_t { return col + first_col; } );
    } else {
        const auto cols_unq = d_cols_unq;
        std::transform( &d_cols_loc[start],
                        &d_cols_loc[end],
                        cols.begin(),
                        [&]( lidx_t col ) -> size_t { return cols_unq[col]; } );
    }

    return cols;
}

} // namespace AMP::LinearAlgebra

#endif
