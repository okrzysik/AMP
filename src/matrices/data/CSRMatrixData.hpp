#ifndef included_AMP_CSRMatrixData_hpp
#define included_AMP_CSRMatrixData_hpp

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/CSRMatrixParameters.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Utilities.h"

#include <memory>
#include <numeric>

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
CSRMatrixData<Policy>::CSRMatrixData( std::shared_ptr<MatrixParametersBase> params )
    : MatrixData( params )
{
    AMPManager::incrementResource( "CSRMatrixData" );
    auto csrParams = std::dynamic_pointer_cast<CSRMatrixParameters<Policy>>( d_pParameters );
    if ( csrParams ) {
        d_is_square   = csrParams->d_is_square;
        d_first_row   = csrParams->d_first_row;
        d_last_row    = csrParams->d_last_row;
        d_first_col   = csrParams->d_first_col;
        d_last_col    = csrParams->d_last_col;
        d_cols        = csrParams->d_cols;
        d_nnz_per_row = csrParams->d_nnz_per_row;
        d_coeffs      = csrParams->d_coeffs;

    } else {
        AMP_ERROR( "Requires CSRParameter object at present" );
    }

    auto memType = AMP::Utilities::getMemoryType( d_cols );
    // the next line should probably not allow for unregistered
    if ( memType == AMP::Utilities::MemoryType::host ||
         memType == AMP::Utilities::MemoryType::unregistered ) {
        size_t N         = d_last_row - d_first_row;
        const size_t nnz = std::accumulate( d_nnz_per_row, d_nnz_per_row + N, 0 );
        std::vector<size_t> remote_dofs;
        for ( auto i = 0u; i < nnz; ++i ) {
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
    // the next line should probably not allow for unregistered
    if ( memType < AMP::Utilities::MemoryType::device ) {
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

        if constexpr ( std::is_same_v<double, scalar_t> ) {
            std::copy( &d_coeffs[offset], &d_coeffs[offset] + n, std::back_inserter( values ) );
        } else {
            std::transform( &d_coeffs[offset],
                            &d_coeffs[offset] + n,
                            std::back_inserter( values ),
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
    // the next line should probably not allow for unregistered
    if ( memType == AMP::Utilities::MemoryType::host ||
         memType == AMP::Utilities::MemoryType::unregistered ) {

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
