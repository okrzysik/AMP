#ifndef included_AMP_CSRMatrixData_hpp
#define included_AMP_CSRMatrixData_hpp

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/CSRMatrixParameters.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Utilities.h"

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
    AMP_ASSERT( AMP::Utilities::getMemoryType( nnz_per_row ) <= AMP::Utilities::MemoryType::host );
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
        // add check for memory location etc and migrate if necessary
        d_is_square     = csrParams->d_is_square;
        d_first_row     = csrParams->d_first_row;
        d_last_row      = csrParams->d_last_row;
        d_first_col     = csrParams->d_first_col;
        d_last_col      = csrParams->d_last_col;
        d_cols          = csrParams->d_cols;
        d_nnz_per_row   = csrParams->d_nnz_per_row;
        d_coeffs        = csrParams->d_coeffs;
        size_t N        = d_last_row - d_first_row;
        d_manage_cols   = false;
        d_manage_nnz    = false;
        d_manage_coeffs = false;

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

        d_leftDOFManager  = matParams->getLeftDOFManager();
        d_rightDOFManager = matParams->getRightDOFManager();
        AMP_ASSERT( d_leftDOFManager && d_rightDOFManager );

        d_is_square = ( d_leftDOFManager->numGlobalDOF() == d_rightDOFManager->numGlobalDOF() );
        d_first_row = d_leftDOFManager->beginDOF();
        d_last_row  = d_leftDOFManager->endDOF();
        d_first_col = d_rightDOFManager->beginDOF();
        d_last_col  = d_rightDOFManager->endDOF();

        auto *nnzPerRow = matParams->entryList();
        auto &cols      = matParams->getColumns();

        if ( d_memory_location <= AMP::Utilities::MemoryType::host ) {

            d_manage_nnz    = false;
            d_manage_coeffs = true;
            d_nnz_per_row   = nnzPerRow;
            d_nnz           = cols.size();

            if constexpr ( std::is_same_v<decltype( d_cols ), decltype( cols.data() )> ) {
                d_manage_cols = false;
                d_cols        = cols.data();
            } else {
                d_manage_cols = true;
                std::allocator<gidx_t> allocator_g;
                d_cols = allocator_g.allocate( d_nnz );
                std::transform(
                    cols.begin(), cols.end(), d_cols, []( size_t col ) -> gidx_t { return col; } );
            }

            d_manage_coeffs = true;
            std::allocator<scalar_t> allocator_s;
            d_coeffs = allocator_s.allocate( d_nnz );

        } else if ( d_memory_location == AMP::Utilities::MemoryType::managed ) {
            AMP_ERROR( "CSRMatrixData: managed memory handling has not been implemented as yet" );
        } else {
            AMP_ERROR( "CSRMatrixData: device memory handling has not been implemented as yet" );
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

    } else {
        AMP_ERROR( "Only host memory deallocation handled at present" );
    }
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
    AMP_ERROR( "Not implemented" );
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
        AMP_ERROR( "CSRMatrixData:getRowByGlobalID not implemented for device memory" );
    }
}

template<typename Policy>
void CSRMatrixData<Policy>::addValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, void *values, const typeID &id )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixData<Policy>::setValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, void *values, const typeID &id )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixData<Policy>::getValuesByGlobalID( size_t num_rows,
                                                 size_t num_cols,
                                                 size_t *rows,
                                                 size_t *cols,
                                                 void *values,
                                                 const typeID &id ) const
{
    AMP_ERROR( "Not implemented" );
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
        AMP_ERROR( "CSRMatrixData:getRowByGlobalID not implemented for device memory" );
    }
}

template<typename Policy>
void CSRMatrixData<Policy>::makeConsistent()
{
    AMP_ERROR( "Not implemented" );
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
