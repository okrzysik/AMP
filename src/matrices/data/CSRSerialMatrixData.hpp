#ifndef included_AMP_CSRSerialMatrixData_hpp
#define included_AMP_CSRSerialMatrixData_hpp

#include "AMP/AMP_TPLs.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/CSRSerialMatrixParameters.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/matrices/data/CSRSerialMatrixData.h"
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
CSRSerialMatrixData<Policy>::CSRSerialMatrixData()
{
    AMPManager::incrementResource( "CSRSerialMatrixData" );
}

template<typename T, template<typename> typename Allocator>
static T *allocate( size_t N )
{
    Allocator<T> alloc;
    return alloc.allocate( N );
}

template<typename Policy>
bool isColValid( typename Policy::gidx_t col, bool is_diag,
		 typename Policy::gidx_t first_row, typename Policy::gidx_t last_row )
{
    bool dValid = is_diag && ( first_row <= col && col < last_row );
    bool odValid = !is_diag && ( col < first_row || last_row <= col);
    return ( dValid || odValid );
}

template<typename Policy>
CSRSerialMatrixData<Policy>::CSRSerialMatrixData( std::shared_ptr<MatrixParametersBase> params,
						  bool is_diag )
{
    AMPManager::incrementResource( "CSRSerialMatrixData" );
    d_pParameters = params;
    auto csrParams = std::dynamic_pointer_cast<CSRSerialMatrixParameters<Policy>>( d_pParameters );
    auto matParams = std ::dynamic_pointer_cast<MatrixParameters>( d_pParameters );

    d_memory_location = d_pParameters->d_memory_location;

    if ( csrParams ) {
        // add check for memory location etc and migrate if necessary
        d_is_diag       = csrParams->d_is_diag;
        d_is_empty      = csrParams->d_is_empty;
        d_first_row     = csrParams->d_first_row;
        d_last_row      = csrParams->d_last_row;
        d_first_col     = csrParams->d_first_col;
        d_last_col      = csrParams->d_last_col;
        d_cols          = csrParams->d_cols;
        d_cols_loc      = csrParams->d_cols_loc;
        d_nnz_per_row   = csrParams->d_nnz_per_row;
        d_coeffs        = csrParams->d_coeffs;
        d_manage_cols   = false;
        d_manage_nnz    = false;
        d_manage_coeffs = false;

        AMP_ERROR( "Should this even be allowed, need to know more about DOFManagers. " );

    } else if ( matParams ) {

        // for now all matrix parameter data is assumed to be on host

        auto leftDOFManager  = matParams->getLeftDOFManager();
        auto rightDOFManager = matParams->getRightDOFManager();
        AMP_ASSERT( leftDOFManager && rightDOFManager );

        d_is_diag   = is_diag;
	d_is_empty  = false;
        d_first_row = leftDOFManager->beginDOF();
        d_last_row  = leftDOFManager->endDOF();
        d_first_col = rightDOFManager->beginDOF();
        d_last_col  = rightDOFManager->endDOF();
	
        size_t nRows       = d_last_row - d_first_row;
        auto *nnzPerRowAll = matParams->entryList();
        auto &cols         = matParams->getColumns();
	AMP_INSIST( !cols.empty(),
		    "CSRSerialMatrixData not constructable from MatrixParameters with emtpy columns" );

	// Count number of nonzeros depeding on value of isDiag
        d_nnz = 0;
	for ( size_t i = 0; i < cols.size(); ++i ) {
 	    if ( isColValid<Policy>( cols[i], d_is_diag, d_first_row, d_last_row ) ) {
	        d_nnz++;
	    }
	}

	// bail out for degenerate case with no nnz
	// may happen in off-diagonal blocks
	if ( d_nnz == 0 ) {
	    d_is_empty = true;
            d_manage_nnz    = false;
            d_manage_coeffs = false;
            d_manage_cols   = false;
	    return;
	}

        if ( d_memory_location <= AMP::Utilities::MemoryType::host ) {

            d_manage_nnz    = true;
            d_manage_coeffs = true;
            d_manage_cols   = true;
	    
            d_nnz_per_row = allocate<lidx_t, std::allocator>( nRows );
            d_cols        = allocate<gidx_t, std::allocator>( d_nnz );
            d_cols_loc    = allocate<lidx_t, std::allocator>( d_nnz );
            d_coeffs      = allocate<scalar_t, std::allocator>( d_nnz );
            d_row_starts  = allocate<lidx_t, std::allocator>( nRows + 1 );

	    // Fill cols and nnz based on local row extents and on/off diag status
	    size_t cgi = 0, cli = 0; // indices into global and local arrays of columns
	    lidx_t nnzFilled = 0;
	    for ( size_t i = 0; i < nRows; ++i ) {
	        d_nnz_per_row[i] = 0;
	        for ( lidx_t j = 0; j < nnzPerRowAll[i]; ++j ) {
		    auto col = cols[cgi++];
		    if ( isColValid<Policy>( col, d_is_diag, d_first_row, d_last_row ) ) {
		        d_nnz_per_row[i]++;
		        d_cols[cli] = col;
			if ( d_is_diag ) {
			    d_cols_loc[cli] = static_cast<lidx_t>(col - d_first_row);
			} else {
			    d_cols_loc[cli] = static_cast<lidx_t>(matParams->d_CommListRight->getLocalGhostID(col));
			}
			cli++;
		        nnzFilled++;
		    }
	        }
	    }
	    AMP_ASSERT( nnzFilled == d_nnz );

        } else if ( ( d_memory_location == AMP::Utilities::MemoryType::managed ) ||
                    ( d_memory_location == AMP::Utilities::MemoryType::device ) ) {
            AMP_ERROR( "CSRSerialMatrixData: managed and device memory support not implemented yet" );
        } else {
            AMP_ERROR( "CSRSerialMatrixData: memory space undefined" );
        }

        if ( d_memory_location < AMP::Utilities::MemoryType::device ) {
            std::exclusive_scan( d_nnz_per_row, d_nnz_per_row + nRows, d_row_starts, 0 );
            d_row_starts[nRows] = d_row_starts[nRows - 1] + d_nnz_per_row[nRows - 1];
        } else {
            AMP_ERROR( "CSRSerialMatrixData: row starts not implemented" );
        }

    } else {
        AMP_ERROR( "Check supplied MatrixParameter object" );
    }
}

template<typename Policy>
CSRSerialMatrixData<Policy>::~CSRSerialMatrixData()
{
    AMPManager::decrementResource( "CSRSerialMatrixData" );
    auto matParams = std ::dynamic_pointer_cast<MatrixParameters>( d_pParameters );

    if ( matParams ) {
        // tackle this case for now
        if ( d_memory_location <= AMP::Utilities::MemoryType::host ) {

            if ( d_row_starts ) {
                std::allocator<lidx_t> allocator_l;
                allocator_l.deallocate( d_row_starts, d_last_row - d_first_row + 1 );
            }
            if ( d_manage_cols ) {
                std::allocator<gidx_t> allocator_g;
                std::allocator<lidx_t> allocator_l;
                allocator_g.deallocate( d_cols, d_nnz );
                allocator_l.deallocate( d_cols_loc, d_nnz );
            }

            if ( d_manage_nnz ) {
                std::allocator<lidx_t> allocator_l;
                allocator_l.deallocate( d_nnz_per_row, d_last_row - d_first_row );
            }

            if ( d_manage_coeffs ) {
                std::allocator<scalar_t> allocator_s;
                allocator_s.deallocate( d_coeffs, d_nnz );
            }
        } else if ( ( d_memory_location == AMP::Utilities::MemoryType::managed ) ||
                    ( d_memory_location == AMP::Utilities::MemoryType::device ) ) {
            AMP_ERROR( "CSRSerialMatrixData: managed and device memory support not implemented yet" );
        } else {
            AMP_ERROR( "CSRSerialMatrixData: memory space undefined" );
        }
    }
}

template<typename Policy>
std::shared_ptr<MatrixData> CSRSerialMatrixData<Policy>::cloneMatrixData() const
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
std::shared_ptr<MatrixData> CSRSerialMatrixData<Policy>::transpose() const
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRSerialMatrixData<Policy>::extractDiagonal( std::shared_ptr<Vector> buf ) const
{
    AMP_ASSERT( buf && buf->numberOfDataBlocks() == 1 ); // temporary constraint
    AMP_ASSERT( buf->isType<scalar_t>( 0 ) );
    AMP_INSIST( d_is_diag,
		"CSRSerialMatrixData<Policy>::extractDiagonal can not be called on off-diag matrices" );

    auto *rawVecData = buf->getRawDataBlock<scalar_t>();
    auto memType     = AMP::Utilities::getMemoryType( rawVecData );
    if ( memType < AMP::Utilities::MemoryType::device ) {

        const size_t N = d_last_row - d_first_row;
        for ( size_t i = 0; i < N; ++i ) {
            const auto start = d_row_starts[i];
            const auto end   = d_row_starts[i + 1];
            // colums are unordered at present
            for ( lidx_t j = start; j < end; ++j ) {
                if ( d_cols[j] == static_cast<gidx_t>( d_first_col + i ) ) {
                    rawVecData[i] = d_coeffs[j];
                    break;
                }
            }
        }
    } else {
        AMP_ERROR( "CSRSerialMatrixData<Policy>::extractDiagonal not implemented for vec and matrix in "
                   "different memory spaces" );
    }
}
template<typename Policy>
void CSRSerialMatrixData<Policy>::getRowByGlobalID( size_t row,
						    std::vector<size_t> &cols,
						    std::vector<double> &values ) const
{
    // Don't do anything on empty matrices
    if ( d_is_empty ) { return; }
    
    AMP_INSIST( row >= static_cast<size_t>( d_first_row ) &&
                    row < static_cast<size_t>( d_last_row ),
                "row must be owned by rank" );
    auto memType = AMP::Utilities::getMemoryType( d_cols );

    if ( memType < AMP::Utilities::MemoryType::device ) {
        const auto row_offset = static_cast<size_t>( row - d_first_row );
        const auto offset     = std::accumulate( d_nnz_per_row, d_nnz_per_row + row_offset, 0 );
        const auto n          = d_nnz_per_row[row_offset];

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

template<typename Policy>
void CSRSerialMatrixData<Policy>::addValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, void *vals, const typeID &id )
{
    // Don't do anything on empty matrices
    if ( d_is_empty ) { return; }
    
    if ( getTypeID<scalar_t>() == id ) {

        if ( d_memory_location < AMP::Utilities::MemoryType::device ) {

            auto values = reinterpret_cast<const scalar_t *>( vals );

            for ( size_t i = 0u; i != num_rows; i++ ) {
                if ( rows[i] >= static_cast<size_t>( d_first_row ) &&
                     rows[i] < static_cast<size_t>( d_last_row ) ) {

                    const auto local_row = rows[i] - d_first_row;
                    const auto start     = d_row_starts[local_row];
                    const auto end       = d_row_starts[local_row + 1];
                    // Inefficient because we don't assume order
                    // not sure it's worth optimizing for our use cases
                    for ( size_t icol = 0; icol < num_cols; ++icol ) {
                        for ( lidx_t j = start; j < end; ++j ) {
                            if ( d_cols[j] == static_cast<gidx_t>( cols[icol] ) ) {
                                d_coeffs[j] += values[num_cols * i + icol];
                            }
                        }
                    }
                } else {
		    AMP_ERROR( "What is d_other_data for?" );
                    // for ( size_t icol = 0; icol < num_cols; ++icol ) {
                    //     d_other_data[rows[i]][cols[icol]] += values[num_cols * i + icol];
                    // }
                }
            }

        } else {
            AMP_ERROR( "CSRSerialMatrixData::addValuesByGlobalID not implemented for device memory" );
        }
    } else {
        AMP_ERROR( "Conversion not implemented" );
    }
}

template<typename Policy>
void CSRSerialMatrixData<Policy>::setValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, void *vals, const typeID &id )
{
    // Don't do anything on empty matrices
    if ( d_is_empty ) { return; }
    
    if ( getTypeID<scalar_t>() == id ) {
        if ( d_memory_location < AMP::Utilities::MemoryType::device ) {

            auto values = reinterpret_cast<const scalar_t *>( vals );

            for ( size_t i = 0u; i != num_rows; i++ ) {

                if ( rows[i] >= static_cast<size_t>( d_first_row ) &&
                     rows[i] < static_cast<size_t>( d_last_row ) ) {
                    const auto local_row = rows[i] - d_first_row;
                    const auto start     = d_row_starts[local_row];
                    const auto end       = d_row_starts[local_row + 1];
		    
                    // Inefficient because we don't assume order
                    // not sure it's worth optimizing for our use cases
                    for ( size_t icol = 0; icol < num_cols; ++icol ) {
                        for ( lidx_t j = start; j < end; ++j ) {
                            if ( d_cols[j] == static_cast<gidx_t>( cols[icol] ) ) {
                                d_coeffs[j] = values[num_cols * i + icol];
                            }
                        }
                    }

                } else {
		    AMP_ERROR( "What is d_ghost_data for?" );
                    // for ( size_t icol = 0; icol < num_cols; ++icol ) {
                    //     d_ghost_data[rows[i]][cols[icol]] = values[num_cols * i + icol];
                    // }
                }
            }

        } else {
            AMP_ERROR( "CSRSerialMatrixData::addValuesByGlobalID not implemented for device memory" );
        }
    } else {
        AMP_ERROR( "Conversion not implemented" );
    }
}

template<typename Policy>
void CSRSerialMatrixData<Policy>::getValuesByGlobalID( size_t num_rows,
                                                 size_t num_cols,
                                                 size_t *rows,
                                                 size_t *cols,
                                                 void *values,
                                                 const typeID &id ) const
{
    // Don't do anything on empty matrices
    if ( d_is_empty ) { return; }
    
    if ( getTypeID<scalar_t>() == id ) {
        if ( d_memory_location < AMP::Utilities::MemoryType::device ) {

            if ( num_rows == 1 && num_cols == 1 ) {

                const auto local_row = rows[0] - d_first_row;
                const auto start     = d_row_starts[local_row];
                const auto end       = d_row_starts[local_row + 1];

                for ( lidx_t i = start; i < end; ++i ) {
                    if ( d_cols[i] == static_cast<gidx_t>( cols[0] ) ) {
                        *( reinterpret_cast<scalar_t *>( values ) ) = d_coeffs[i];
                    }
                }
            } else {
                AMP_ERROR( "CSRSerialMatrixData::getValuesByGlobalID not implemented for num_rows>1 || "
                           "num_cols > 1" );
            }

        } else {
            AMP_ERROR( "CSRSerialMatrixData::getValuesByGlobalID not implemented for device memory" );
        }
    } else {
        AMP_ERROR( "Not implemented" );
    }
}

template<typename Policy>
std::vector<size_t> CSRSerialMatrixData<Policy>::getColumnIDs( size_t row ) const
{
    // Don't do anything on empty matrices
  if ( d_is_empty ) { return std::vector<size_t>(); }
    
    AMP_INSIST( row >= static_cast<size_t>( d_first_row ) &&
                    row < static_cast<size_t>( d_last_row ),
                "row must be owned by rank" );
    AMP_INSIST( d_cols && d_nnz_per_row, "Must be initialized" );

    
    auto memType = AMP::Utilities::getMemoryType( d_cols );

    if ( memType < AMP::Utilities::MemoryType::device ) {

        std::vector<size_t> cols;
        const auto row_offset = static_cast<size_t>( row - d_first_row );
        const auto offset     = std::accumulate( d_nnz_per_row, d_nnz_per_row + row_offset, 0 );
        const auto n          = d_nnz_per_row[row_offset];

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

/********************************************************
 * Get the number of rows/columns in the matrix          *
 ********************************************************/
template<typename Policy>
size_t CSRSerialMatrixData<Policy>::numLocalRows() const
{
    return static_cast<size_t>( d_last_row - d_first_row );
}

template<typename Policy>
size_t CSRSerialMatrixData<Policy>::numLocalColumns() const
{
    return static_cast<size_t>( d_last_col - d_first_col );
}

/********************************************************
 * Get iterators                                         *
 ********************************************************/
template<typename Policy>
size_t CSRSerialMatrixData<Policy>::beginRow() const
{
    return static_cast<size_t>( d_first_row );
}

template<typename Policy>
size_t CSRSerialMatrixData<Policy>::endRow() const
{
    return static_cast<size_t>( d_last_row );
}

} // namespace AMP::LinearAlgebra

#endif
