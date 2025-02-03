#ifndef included_AMP_CSRLocalMatrixData_hpp
#define included_AMP_CSRLocalMatrixData_hpp

#include "AMP/AMP_TPLs.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/CSRMatrixParameters.h"
#include "AMP/matrices/MatrixParameters.h"
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
    d_pParameters  = params;
    auto csrParams = std::dynamic_pointer_cast<CSRMatrixParameters<Policy>>( d_pParameters );
    auto matParams = std ::dynamic_pointer_cast<MatrixParameters>( d_pParameters );

    if ( csrParams ) {
        // Pull out block specific parameters
        auto &blParams = d_is_diag ? csrParams->d_diag : csrParams->d_off_diag;

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

        const auto &getRow = matParams->getRowFunction();

        if ( !getRow ) {
            // Without a getRow function there is no way to set the nnz structure
            // This route will just create an empty matrix
            // This can happen e.g. when constructing the output matrix of SpGEMM
            d_nnz      = 0;
            d_is_empty = true;
            return;
        }

        // Count number of nonzeros
        d_nnz = 0;
        for ( gidx_t i = d_first_row; i < d_last_row; ++i ) {
            for ( auto &&col : getRow( i ) ) {
                if ( isColValid<Policy>( col, d_is_diag, d_first_col, d_last_col ) ) {
                    ++d_nnz;
                }
            }
        }

        // bail out for degenerate case with no nnz
        // may happen in off-diagonal blocks
        if ( d_nnz == 0 ) {
            d_is_empty = true;
            return;
        }
        d_is_empty = false;

        // Allocate internal arrays
        d_row_starts = sharedArrayBuilder( d_num_rows + 1, d_lidxAllocator );
        d_cols       = sharedArrayBuilder( d_nnz, d_gidxAllocator );
        d_coeffs     = sharedArrayBuilder( d_nnz, d_scalarAllocator );

        // Fill cols and nnz based on local row extents and on/off diag status
        lidx_t nnzFilled = 0;
        d_row_starts[0]  = 0;
        for ( lidx_t row = 0; row < d_num_rows; ++row ) {
            auto cols = getRow( d_first_row + row );
            // initialize next rs to this one and push forward as nz's added to this row
            d_row_starts[row + 1] = d_row_starts[row];
            for ( auto &&col : cols ) {
                if ( isColValid<Policy>( col, d_is_diag, d_first_col, d_last_col ) ) {
                    d_row_starts[row + 1]++;
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
        // an encasing CSRMatrixData object). This is used for remote blocks
        // in SpGEMM
        d_nnz      = 0;
        d_is_empty = true;
        return;
    }

    // now do setup that is independent of MatrixParameter datatype
    d_max_row_len = 0;
    for ( lidx_t row = 0; row < d_num_rows; ++row ) {
        const auto ncols = d_row_starts[row + 1] - d_row_starts[row];
        d_max_row_len    = d_max_row_len < ncols ? ncols : d_max_row_len;
    }

    // local columns always owned internally
    d_cols_loc = sharedArrayBuilder( d_nnz, d_lidxAllocator );

    // fill in local column indices
    globalToLocalColumns();
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::sortColumns( MatrixSortScheme sort_type )
{
    typedef std::tuple<lidx_t, gidx_t, scalar_t> tuple_t;

    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRSerialMatrixData::sortColumns not implemented for device memory" );
    AMP_INSIST( sort_type <= MatrixSortScheme::hypre, "Unknown sorting type requested" );

    if ( d_is_empty || sort_type == MatrixSortScheme::unsorted ) {
        return;
    }

    std::vector<tuple_t> rTpl( d_max_row_len );
    for ( lidx_t row = 0; row < d_num_rows; ++row ) {
        const auto rs      = d_row_starts[row];
        const auto row_len = d_row_starts[row + 1] - rs;

        for ( lidx_t k = 0; k < row_len; ++k ) {
            rTpl[k] = std::make_tuple( d_cols_loc[rs + k], d_cols[rs + k], d_coeffs[rs + k] );
        }

        if ( sort_type == MatrixSortScheme::hypre && d_is_diag ) {
            std::sort( rTpl.begin(),
                       rTpl.begin() + row_len,
                       [row]( const tuple_t &a, const tuple_t &b ) -> bool {
                           const lidx_t lca = std::get<0>( a ), lcb = std::get<0>( b );
                           return row != lcb && ( lca < lcb || lca == row );
                       } );
        } else if ( sort_type == MatrixSortScheme::ascending ||
                    sort_type == MatrixSortScheme::hypre ) {
            // note: ascending sort also for offd hypre convention
            std::sort( rTpl.begin(),
                       rTpl.begin() + row_len,
                       []( const tuple_t &a, const tuple_t &b ) -> bool {
                           // note that this sorts on global index, not local
                           return std::get<1>( a ) < std::get<1>( b );
                       } );
        }

        for ( lidx_t k = 0; k < row_len; ++k ) {
            d_cols_loc[rs + k] = std::get<0>( rTpl[k] );
            d_cols[rs + k]     = std::get<1>( rTpl[k] );
            d_coeffs[rs + k]   = std::get<2>( rTpl[k] );
        }
    }

    // re-write column map to match now permuted columns
    if ( !d_is_diag && sort_type != MatrixSortScheme::hypre ) {
        for ( lidx_t n = 0; n < d_nnz; ++n ) {
            d_cols_unq[d_cols_loc[n]] = d_cols[n];
        }
    } else if ( !d_is_diag ) {
        // hypre requires sorted column map, so sort it and
        // instead change local indices to match
        std::sort( d_cols_unq.get(), d_cols_unq.get() + d_ncols_unq );
        for ( lidx_t n = 0; n < d_nnz; ++n ) {
            auto it =
                std::lower_bound( d_cols_unq.get(), d_cols_unq.get() + d_ncols_unq, d_cols[n] );
            d_cols_loc[n] = static_cast<lidx_t>( std::distance( d_cols_unq.get(), it ) );
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

    cloneData = std::make_shared<CSRLocalMatrixData>( d_pParameters,
                                                      d_memory_location,
                                                      d_first_row,
                                                      d_last_row,
                                                      d_first_col,
                                                      d_last_col,
                                                      d_is_diag );

    cloneData->d_is_empty = d_is_empty;
    cloneData->d_nnz      = d_nnz;

    if ( !d_is_empty ) {
        cloneData->d_row_starts = sharedArrayBuilder( d_num_rows + 1, d_lidxAllocator );
        cloneData->d_cols       = sharedArrayBuilder( d_nnz, d_gidxAllocator );
        cloneData->d_cols_loc   = sharedArrayBuilder( d_nnz, d_lidxAllocator );
        cloneData->d_coeffs     = sharedArrayBuilder( d_nnz, d_scalarAllocator );

        AMP::Utilities::Algorithms<lidx_t>::copy_n(
            d_row_starts.get(), d_num_rows + 1, cloneData->d_row_starts.get() );
        AMP::Utilities::Algorithms<gidx_t>::copy_n( d_cols.get(), d_nnz, cloneData->d_cols.get() );
        AMP::Utilities::Algorithms<lidx_t>::copy_n(
            d_cols_loc.get(), d_nnz, cloneData->d_cols_loc.get() );
        AMP::Utilities::Algorithms<scalar_t>::fill_n( cloneData->d_coeffs.get(), d_nnz, 0.0 );

    } else {
        cloneData->d_row_starts = nullptr;
        cloneData->d_cols       = nullptr;
        cloneData->d_cols_loc   = nullptr;
        cloneData->d_coeffs     = nullptr;
    }

    return cloneData;
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::setNNZ( const std::vector<lidx_t> &nnz )
{
    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRLocalMatrixData::setNNZ not implemented on device yet" );

    // allocate and fill rowstarts from scan of passed nnz vector
    d_row_starts = sharedArrayBuilder( d_num_rows + 1, d_lidxAllocator );
    std::exclusive_scan( nnz.begin(), nnz.end(), d_row_starts.get(), 0 );
    d_row_starts[d_num_rows] = d_row_starts[d_num_rows - 1] + nnz[d_num_rows - 1];

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
    std::fill( d_coeffs.get(), d_coeffs.get() + d_nnz, 0.0 );
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::globalToLocalColumns()
{
    if ( d_is_empty ) {
        return;
    }

    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRLocalMatrixData::globalToLocalColumns not implemented on device yet" );

    if ( d_is_diag ) {
        for ( lidx_t n = 0; n < d_nnz; ++n ) {
            d_cols_loc[n] = static_cast<lidx_t>( d_cols[n] - d_first_col );
        }
    } else {
        // for offd setup column map as part of the process
        std::unordered_map<gidx_t, lidx_t> colMap;
        d_ncols_unq = 0;
        for ( lidx_t n = 0; n < d_nnz; ++n ) {
            const auto inserted = colMap.insert( { d_cols[n], d_ncols_unq } );
            if ( inserted.second ) {
                ++d_ncols_unq;
            }
            d_cols_loc[n] = inserted.first->second;
        }

        // offd column map also always owned internally
        d_cols_unq = sharedArrayBuilder( d_ncols_unq, d_gidxAllocator );
        for ( auto it = colMap.begin(); it != colMap.end(); ++it ) {
            d_cols_unq[it->second] = it->first;
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

    const auto offset = d_row_starts[local_row];
    auto n            = d_row_starts[local_row + 1] - d_row_starts[local_row];

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
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::getValuesByGlobalID( const size_t local_row,
                                                                 const size_t col,
                                                                 void *values,
                                                                 const typeID &id ) const
{
    // Don't do anything on empty matrices
    if ( d_is_empty ) {
        return;
    }

    AMP_INSIST( getTypeID<scalar_t>() == id,
                "CSRLocalMatrixData::getValuesByGlobalID called with inconsistent typeID" );

    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRLocalMatrixData::getValuesByGlobalID not implemented for device memory" );

    const auto start = d_row_starts[local_row];
    auto end         = d_row_starts[local_row + 1];

    for ( lidx_t i = start; i < end; ++i ) {
        if ( d_cols[i] == static_cast<gidx_t>( col ) ) {
            *( reinterpret_cast<scalar_t *>( values ) ) = d_coeffs[i];
        }
    }
}

template<typename Policy, class Allocator>
void CSRLocalMatrixData<Policy, Allocator>::addValuesByGlobalID( const size_t num_cols,
                                                                 const size_t local_row,
                                                                 const size_t *cols,
                                                                 const scalar_t *vals,
                                                                 const typeID &id )
{
    if ( d_is_empty ) {
        return;
    }

    AMP_INSIST( getTypeID<scalar_t>() == id,
                "CSRLocalMatrixData::addValuesByGlobalID called with inconsistent typeID" );

    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRLocalMatrixData::addValuesByGlobalID not implemented for device memory" );

    const auto start = d_row_starts[local_row];
    auto end         = d_row_starts[local_row + 1];

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
void CSRLocalMatrixData<Policy, Allocator>::setValuesByGlobalID( const size_t num_cols,
                                                                 const size_t local_row,
                                                                 const size_t *cols,
                                                                 const scalar_t *vals,
                                                                 const typeID &id )
{
    if ( d_is_empty ) {
        return;
    }

    AMP_INSIST( getTypeID<scalar_t>() == id,
                "CSRLocalMatrixData::setValuesByGlobalID called with inconsistent typeID" );

    AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                "CSRLocalMatrixData::setValuesByGlobalID not implemented for device memory" );

    const auto start = d_row_starts[local_row];
    auto end         = d_row_starts[local_row + 1];

    // Inefficient because we don't assume order
    // not sure it's worth optimizing for our use cases
    for ( size_t icol = 0; icol < num_cols; ++icol ) {
        for ( lidx_t j = start; j < end; ++j ) {
            if ( d_cols[j] == static_cast<gidx_t>( cols[icol] ) ) {
                d_coeffs[j] = vals[icol];
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

    AMP_INSIST( d_cols && d_row_starts,
                "CSRLocalMatrixData::getColumnIDs nnz layout must be initialized" );

    std::vector<size_t> cols;
    const auto offset = d_row_starts[local_row];
    auto n            = d_row_starts[local_row + 1] - d_row_starts[local_row];

    if constexpr ( std::is_same_v<size_t, gidx_t> ) {
        std::copy( &d_cols[offset], &d_cols[offset] + n, std::back_inserter( cols ) );
    } else {
        std::transform( &d_cols[offset],
                        &d_cols[offset] + n,
                        std::back_inserter( cols ),
                        []( size_t col ) -> gidx_t { return col; } );
    }
    return cols;
}

} // namespace AMP::LinearAlgebra

#endif
