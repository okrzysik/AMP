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

template<typename T, template<typename> typename Allocator>
static T *allocate( size_t N )
{
    Allocator<T> alloc;
    return alloc.allocate( N );
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
        // add check for memory location etc and migrate if necessary
        d_is_square     = csrParams->d_is_square;
        d_first_row     = csrParams->d_first_row;
        d_last_row      = csrParams->d_last_row;
        d_first_col     = csrParams->d_first_col;
        d_last_col      = csrParams->d_last_col;
        d_cols          = csrParams->d_cols;
        d_nnz_per_row   = csrParams->d_nnz_per_row;
        d_coeffs        = csrParams->d_coeffs;
        d_manage_cols   = false;
        d_manage_nnz    = false;
        d_manage_coeffs = false;

        size_t N = d_last_row - d_first_row;

        auto memType = AMP::Utilities::getMemoryType( d_cols );

        if ( memType != AMP::Utilities::MemoryType::device ) {
            d_nnz = NNZ<Policy>( N, d_nnz_per_row );
            std::vector<size_t> remote_dofs;
            for ( lidx_t i = 0; i < d_nnz; ++i ) {
                if ( ( d_cols[i] < d_first_col ) || ( d_cols[i] >= d_last_col ) ) {
                    remote_dofs.push_back( d_cols[i] );
                }
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
            AMP_WARNING( "CSRMatrixData: device memory handling has not been implemented as yet" );
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
        size_t N    = d_last_row - d_first_row;

        auto *nnzPerRow = matParams->entryList();
        auto &cols      = matParams->getColumns();
        d_nnz           = cols.empty() ? NNZ<Policy>( N, nnzPerRow ) : cols.size();

        if ( d_memory_location <= AMP::Utilities::MemoryType::host ) {

            d_manage_nnz    = false;
            d_manage_coeffs = true;
            d_manage_cols   = cols.empty() ? true : false;
            d_nnz_per_row   = nnzPerRow;

            if constexpr ( std::is_same_v<decltype( d_cols ), decltype( cols.data() )> ) {
                if ( d_manage_cols ) {
                    d_cols = allocate<gidx_t, std::allocator>( d_nnz );
                } else {
                    d_cols = cols.data();
                }
            } else {
                d_manage_cols = true;
                d_cols        = allocate<gidx_t, std::allocator>( d_nnz );
                std::transform(
                    cols.begin(), cols.end(), d_cols, []( size_t col ) -> gidx_t { return col; } );
            }

            d_manage_coeffs = true;
            d_coeffs        = allocate<scalar_t, std::allocator>( d_nnz );

            d_row_starts = allocate<lidx_t, std::allocator>( N + 1 );

        } else if ( ( d_memory_location == AMP::Utilities::MemoryType::managed ) ||
                    ( d_memory_location == AMP::Utilities::MemoryType::device ) ) {

            d_manage_cols   = true;
            d_manage_nnz    = true;
            d_manage_coeffs = true;

#ifdef AMP_USE_UMPIRE
            auto &resourceManager = umpire::ResourceManager::getInstance();
            auto allocator        = ( d_memory_location == AMP::Utilities::MemoryType::managed ) ?
                                        resourceManager.getAllocator( "UM" ) :
                                        resourceManager.getAllocator( "DEVICE" );

            d_nnz_per_row = static_cast<lidx_t *>( allocator.allocate( N * sizeof( lidx_t ) ) );
            d_row_starts =
                static_cast<lidx_t *>( allocator.allocate( ( N + 1 ) * sizeof( lidx_t ) ) );

            d_cols   = static_cast<gidx_t *>( allocator.allocate( d_nnz * sizeof( gidx_t ) ) );
            d_coeffs = static_cast<scalar_t *>( allocator.allocate( d_nnz * sizeof( scalar_t ) ) );

            AMP_ASSERT( d_nnz_per_row && d_cols && d_coeffs );

            static_assert( std::is_same_v<decltype( d_nnz_per_row ), decltype( nnzPerRow )> );
            resourceManager.copy( d_nnz_per_row, nnzPerRow );

            if constexpr ( std::is_same_v<decltype( d_cols ), decltype( cols.data() )> ) {
                resourceManager.copy( d_cols, cols.data() );
            } else {
                if ( d_memory_location == AMP::Utilities::MemoryType::managed ) {
                    std::transform( cols.begin(), cols.end(), d_cols, []( size_t col ) -> gidx_t {
                        return col;
                    } );
                } else {
                    AMP_ERROR( "Not implemented" );
                }
            }
#else
            AMP_ERROR(
                "CSRMatrixData: managed and device memory handling without Umpire has not been "
                "implemented as yet" );
#endif
        } else {
            AMP_ERROR( "CSRMatrixData: memory space undefined" );
        }

        if ( d_memory_location < AMP::Utilities::MemoryType::device ) {
            std::exclusive_scan( d_nnz_per_row, d_nnz_per_row + N, d_row_starts, 0 );
            d_row_starts[N] = d_row_starts[N - 1] + d_nnz_per_row[N - 1];
        } else {
            AMP_ERROR( "CSRMatrixData: row starts not implemented" );
        }

    } else {
        AMP_ERROR( "Check supplied MatrixParameter object" );
    }
}

template<typename Policy>
CSRMatrixData<Policy>::~CSRMatrixData()
{
    AMPManager::decrementResource( "CSRMatrixData" );

    // tackle this case for now
    if ( d_memory_location <= AMP::Utilities::MemoryType::host ) {

        if ( d_row_starts ) {
            std::allocator<lidx_t> allocator_l;
            allocator_l.deallocate( d_row_starts, d_last_row - d_first_row + 1 );
        }
        if ( d_manage_cols ) {
            std::allocator<gidx_t> allocator_g;
            allocator_g.deallocate( d_cols, d_nnz );
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

#ifdef AMP_USE_UMPIRE
        auto &resourceManager = umpire::ResourceManager::getInstance();
        auto allocator        = ( d_memory_location == AMP::Utilities::MemoryType::managed ) ?
                                    resourceManager.getAllocator( "UM" ) :
                                    resourceManager.getAllocator( "DEVICE" );

        allocator.deallocate( d_row_starts );

        if ( d_manage_cols )
            allocator.deallocate( d_cols );

        if ( d_manage_nnz )
            allocator.deallocate( d_nnz_per_row );

        if ( d_manage_coeffs )
            allocator.deallocate( d_coeffs );
#else
        AMP_ERROR( "CSRMatrixData: managed and device memory handling without Umpire has not been "
                   "implemented as yet" );
#endif

    } else {
        AMP_ERROR( "CSRMatrixData: memory space undefined" );
    }
}

template<typename Policy>
std::shared_ptr<MatrixData> CSRMatrixData<Policy>::cloneMatrixData() const
{
    std::shared_ptr<CSRMatrixData> cloneData;
#ifdef AMP_USE_UMPIRE
    auto &resourceManager = umpire::ResourceManager::getInstance();

    umpire::Allocator allocator;
    if ( d_memory_location <= AMP::Utilities::MemoryType::host )
        allocator = resourceManager.getAllocator( "HOST" );
    else if ( d_memory_location == AMP::Utilities::MemoryType::managed )
        allocator = resourceManager.getAllocator( "UM" );
    else if ( d_memory_location == AMP::Utilities::MemoryType::device )
        allocator = resourceManager.getAllocator( "DEVICE" );
    else
        AMP_ERROR( "Unsupported memory location" );

    cloneData = std::make_shared<CSRMatrixData<Policy>>();

    cloneData->d_memory_location = d_memory_location;

    cloneData->d_is_square = d_is_square;
    cloneData->d_first_row = d_first_row;
    cloneData->d_last_row  = d_last_row;
    cloneData->d_first_col = d_first_col;
    cloneData->d_last_col  = d_last_col;
    cloneData->d_nnz       = d_nnz;

    cloneData->d_manage_nnz    = true;
    cloneData->d_manage_coeffs = true;
    cloneData->d_manage_cols   = true;

    size_t N = d_last_row - d_first_row;

    cloneData->d_nnz_per_row = static_cast<lidx_t *>( allocator.allocate( N * sizeof( lidx_t ) ) );
    cloneData->d_row_starts =
        static_cast<lidx_t *>( allocator.allocate( ( N + 1 ) * sizeof( lidx_t ) ) );

    cloneData->d_cols = static_cast<gidx_t *>( allocator.allocate( d_nnz * sizeof( gidx_t ) ) );
    cloneData->d_coeffs =
        static_cast<scalar_t *>( allocator.allocate( d_nnz * sizeof( scalar_t ) ) );

    AMP_ERROR( "Incomplete -- add in DOF manager and base members" );
#else
    AMP_ERROR( "CSRMatrixData: managed and device memory handling without Umpire has not been "
               "implemented as yet" );
#endif
    return cloneData;
}

template<typename Policy>
std::shared_ptr<MatrixData> CSRMatrixData<Policy>::transpose() const
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixData<Policy>::extractDiagonal( std::shared_ptr<Vector> buf ) const
{
    AMP_ASSERT( buf && buf->numberOfDataBlocks() == 1 ); // temporary constraint
    AMP_ASSERT( buf->isType<scalar_t>( 0 ) );

    auto *rawVecData  = buf->getRawDataBlock<scalar_t>();
    auto vDataMemType = AMP::Utilities::getMemoryType( rawVecData );
    if ( vDataMemType == d_memory_location ) {
        if ( d_memory_location < AMP::Utilities::MemoryType::device ) {

            const size_t N = d_last_row - d_first_row;
            for ( size_t i = 0; i < N; ++i ) {
                rawVecData[i] = d_coeffs[d_row_starts[i]];
            }

        } else {
            AMP_ERROR(
                "CSRMatrixData<Policy>::extractDiagonal not implemented for vec and matrix in "
                "different memory spaces" );
        }
    }
}
template<typename Policy>
void CSRMatrixData<Policy>::getRowByGlobalID( size_t row,
                                              std::vector<size_t> &cols,
                                              std::vector<double> &values ) const
{
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
        AMP_ERROR( "CSRMatrixData::getRowByGlobalID not implemented for device memory" );
    }
}

template<typename Policy>
void CSRMatrixData<Policy>::addValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, void *values, const typeID &id )
{
    if ( getTypeID<scalar_t>() == id ) {
        if ( d_memory_location < AMP::Utilities::MemoryType::device ) {

            if ( num_rows == 1 ) {

                const auto local_row = ( *rows ) - d_first_row;
                const auto start     = d_row_starts[local_row];
                const auto end       = d_row_starts[local_row + 1];
                // O(nnz_per_row)^2 because we don't assume order
                // not sure it's worth optimizing for our use cases
                for ( size_t icol = 0; icol < num_cols; ++icol ) {
                    for ( lidx_t i = start; i < end; ++i ) {
                        if ( d_cols[i] == static_cast<lidx_t>( cols[icol] ) ) {
                            d_coeffs[i] += ( reinterpret_cast<scalar_t *>( values )[icol] );
                        }
                    }
                }

            } else {
                AMP_ERROR( "CSRMatrixData::addValuesByGlobalID not implemented for num_rows>1" );
            }

        } else {
            AMP_ERROR( "CSRMatrixData::addValuesByGlobalID not implemented for device memory" );
        }
    } else {
        AMP_ERROR( "Not implemented" );
    }
}

template<typename Policy>
void CSRMatrixData<Policy>::setValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, void *values, const typeID &id )
{
    if ( getTypeID<scalar_t>() == id ) {
        if ( d_memory_location < AMP::Utilities::MemoryType::device ) {

            if ( num_rows == 1 ) {

                const auto local_row = ( *rows ) - d_first_row;
                const auto start     = d_row_starts[local_row];
                const auto end       = d_row_starts[local_row + 1];

                // O(nnz_per_row)^2 because we don't assume order
                // not sure it's worth optimizing for our use cases
                for ( size_t icol = 0; icol < num_cols; ++icol ) {
                    for ( lidx_t i = start; i < end; ++i ) {
                        if ( d_cols[i] == static_cast<lidx_t>( cols[icol] ) ) {
                            d_coeffs[i] = ( reinterpret_cast<scalar_t *>( values )[icol] );
                        }
                    }
                }
            } else {
                AMP_ERROR( "CSRMatrixData::addValuesByGlobalID not implemented for num_rows>1 " );
            }

        } else {
            AMP_ERROR( "CSRMatrixData::addValuesByGlobalID not implemented for device memory" );
        }
    } else {
        AMP_ERROR( "Not implemented" );
    }
}

template<typename Policy>
void CSRMatrixData<Policy>::getValuesByGlobalID( size_t num_rows,
                                                 size_t num_cols,
                                                 size_t *rows,
                                                 size_t *cols,
                                                 void *values,
                                                 const typeID &id ) const
{
    if ( getTypeID<scalar_t>() == id ) {
        if ( d_memory_location < AMP::Utilities::MemoryType::device ) {

            if ( num_rows == 1 && num_cols == 1 ) {

                const auto local_row = ( *rows ) - d_first_row;
                const auto start     = d_row_starts[local_row];
                const auto end       = d_row_starts[local_row + 1];

                for ( lidx_t i = start; i < end; ++i ) {
                    if ( d_cols[i] == static_cast<lidx_t>( *cols ) ) {
                        *( reinterpret_cast<scalar_t *>( values ) ) = d_coeffs[i];
                    }
                }
            } else {
                AMP_ERROR( "CSRMatrixData::getValuesByGlobalID not implemented for num_rows>1 || "
                           "num_cols > 1" );
            }

        } else {
            AMP_ERROR( "CSRMatrixData::getValuesByGlobalID not implemented for device memory" );
        }
    } else {
        AMP_ERROR( "Not implemented" );
    }
}

template<typename Policy>
void CSRMatrixData<Policy>::makeConsistent( AMP::LinearAlgebra::ScatterType )
{
    // at present do nothing since the way AMP is dealing with face/vertex values
    // seems to preclude the need at present
}

template<typename Policy>
std::vector<size_t> CSRMatrixData<Policy>::getColumnIDs( size_t row ) const
{
    AMP_INSIST( row >= static_cast<size_t>( d_first_row ) &&
                    row < static_cast<size_t>( d_last_row ),
                "row must be owned by rank" );
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
        AMP_ERROR( "CSRMatrixData:getColumnIDs not implemented for device memory" );
    }
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
