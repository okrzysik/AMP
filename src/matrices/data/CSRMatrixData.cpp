#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/CSRMatrixParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Utilities.h"

#include <memory>
#include <numeric>

namespace AMP::LinearAlgebra {


/********************************************************
 * Constructors/Destructor                               *
 ********************************************************/
CSRMatrixData::CSRMatrixData() { AMPManager::incrementResource( "CSRMatrixData" ); }
CSRMatrixData::CSRMatrixData( std::shared_ptr<MatrixParametersBase> params ) : MatrixData( params )
{
    AMPManager::incrementResource( "CSRMatrixData" );
    auto csrParams = std::dynamic_pointer_cast<CSRMatrixParameters>( d_pParameters );
    if ( csrParams ) {

        d_first_row   = csrParams->d_first_row;
        d_last_row    = csrParams->d_last_row;
        d_cols        = csrParams->d_cols;
        d_nnz_per_row = csrParams->d_nnz_per_row;
        d_cols        = csrParams->d_cols;
        d_coeffs      = csrParams->d_coeffs;

    } else {
        AMP_ERROR( "Requires CSRParameter object at present" );
    }

    if ( AMP::Utilities::getMemoryType( d_cols ) == AMP::Utilities::MemoryType::host ) {
        size_t N         = d_last_row - d_first_row + 1;
        const size_t nnz = std::accumulate( d_nnz_per_row, d_nnz_per_row + N, 0 );
        std::vector<size_t> remote_dofs;
        for ( auto i = 0u; i < nnz; ++i ) {
            if ( ( d_cols[i] < d_first_row ) || ( d_cols[i] > d_last_row ) ) {
                remote_dofs.push_back( d_cols[i] );
            }
        }
        AMP::Utilities::unique( remote_dofs );
        const auto &comm = getComm();
        d_rightDOFManager =
            std::make_shared<AMP::Discretization::DOFManager>( N, comm, remote_dofs );

    } else {
        AMP_ERROR( "device memory handling has not been implemented as yet" );
    }
}

CSRMatrixData::~CSRMatrixData() { AMPManager::decrementResource( "CSRMatrixData" ); }

std::shared_ptr<MatrixData> CSRMatrixData::cloneMatrixData() const
{
    AMP_ERROR( "Not implemented" );
}

std::shared_ptr<MatrixData> CSRMatrixData::transpose() const { AMP_ERROR( "Not implemented" ); }

void CSRMatrixData::extractDiagonal( std::shared_ptr<Vector> buf ) const
{
    AMP_ERROR( "Not implemented" );
}

void CSRMatrixData::getRowByGlobalID( size_t row,
                                      std::vector<size_t> &cols,
                                      std::vector<double> &values ) const
{
    AMP_ERROR( "Not implemented" );
}

void CSRMatrixData::addValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, void *values, const typeID &id )
{
    AMP_ERROR( "Not implemented" );
}

void CSRMatrixData::setValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, void *values, const typeID &id )
{
    AMP_ERROR( "Not implemented" );
}

void CSRMatrixData::getValuesByGlobalID( size_t num_rows,
                                         size_t num_cols,
                                         size_t *rows,
                                         size_t *cols,
                                         void *values,
                                         const typeID &id ) const
{
    AMP_ERROR( "Not implemented" );
}

std::vector<size_t> CSRMatrixData::getColumnIDs( size_t row ) const
{
    AMP_ERROR( "Not implemented" );
}

void CSRMatrixData::makeConsistent() { AMP_ERROR( "Not implemented" ); }

std::shared_ptr<Discretization::DOFManager> CSRMatrixData::getRightDOFManager() const
{
    return d_rightDOFManager;
}

std::shared_ptr<Discretization::DOFManager> CSRMatrixData::getLeftDOFManager() const
{
    return d_leftDOFManager;
}

/********************************************************
 * Get the number of rows/columns in the matrix          *
 ********************************************************/
size_t CSRMatrixData::numLocalRows() const { return d_last_row - d_first_row + 1; }

size_t CSRMatrixData::numGlobalRows() const { AMP_ERROR( "Not implemented" ); }

size_t CSRMatrixData::numLocalColumns() const
{
    AMP_ERROR( "Not implemented" );
    return 0;
}

size_t CSRMatrixData::numGlobalColumns() const { AMP_ERROR( "Not implemented" ); }


/********************************************************
 * Get iterators                                         *
 ********************************************************/
size_t CSRMatrixData::beginRow() const { return d_first_row; }
size_t CSRMatrixData::endRow() const { return d_last_row + 1; }


} // namespace AMP::LinearAlgebra
