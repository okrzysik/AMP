#ifndef included_MatrixDataTransforms_h
#define included_MatrixDataTransforms_h

#include "AMP/matrices/Matrix.h"

#include <numeric>
#include <type_traits>

namespace AMP::LinearAlgebra {

template<typename Policy>
void transformDofToCSR( std::shared_ptr<Matrix> matrix,
                        typename Policy::gidx_t &firstRow,
                        typename Policy::gidx_t &endRow,
                        std::vector<typename Policy::lidx_t> &nnz_d,
                        std::vector<typename Policy::lidx_t> &rowstart_d,
                        std::vector<typename Policy::gidx_t> &cols_d,
                        std::vector<typename Policy::lidx_t> &cols_loc_d,
                        std::vector<typename Policy::scalar_t> &coeffs_d,
                        std::vector<typename Policy::lidx_t> &nnz_od,
                        std::vector<typename Policy::lidx_t> &rowstart_od,
                        std::vector<typename Policy::gidx_t> &cols_od,
                        std::vector<typename Policy::lidx_t> &cols_loc_od,
                        std::vector<typename Policy::scalar_t> &coeffs_od,
			typename Policy::lidx_t &nnz_pad)
{
    using gidx_t   = typename Policy::gidx_t;
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;
  
    AMP_ASSERT( matrix );

    firstRow = static_cast<gidx_t>( matrix->beginRow() );
    endRow   = static_cast<gidx_t>( matrix->endRow() );

    // make a commlist from dof manager in input matrix
    // ripped from VectorBuilder.cpp
    // there must be a better way...
    auto rightDOFs                               = matrix->getRightDOFManager();
    auto remote_DOFs                             = rightDOFs->getRemoteDOFs();
    auto comm                                    = rightDOFs->getComm();
    bool have_ghosts                             = comm.anyReduce( !remote_DOFs.empty() );
    std::shared_ptr<CommunicationList> comm_list = nullptr;
    if ( have_ghosts ) {
        // Construct the communication list
        auto params           = std::make_shared<CommunicationListParameters>();
        params->d_comm        = comm;
        params->d_localsize   = rightDOFs->numLocalDOF();
        params->d_remote_DOFs = remote_DOFs;
        comm_list             = std::make_shared<CommunicationList>( params );
    }
    comm.barrier();

    // loop over rows and examine structure
    // nearly all of this comes straight from the CSRSerialMatrixData constructor...
    for ( auto row = firstRow; row < endRow; ++row ) {

        std::vector<size_t> rcols;
        std::vector<scalar_t> rvals;

        matrix->getRowByGlobalID( row, rcols, rvals );

        // Loop over columns and insert into on and off diagonal blocks
        lidx_t nnzd = 0, nnzod = 0;
        for ( size_t i = 0; i < rcols.size(); ++i ) {
            const auto col = static_cast<gidx_t>( rcols[i] );
            if ( firstRow <= col && col < endRow ) {
                nnzd++;
                cols_d.push_back( col );
                cols_loc_d.push_back( static_cast<lidx_t>( col - firstRow ) );
                coeffs_d.push_back( static_cast<scalar_t>( rvals[i] ) );
            } else if ( have_ghosts ) {
                nnzod++;
                cols_od.push_back( col );
                auto loc_col = comm_list->getLocalGhostID( rcols[i] );
                cols_loc_od.push_back( static_cast<lidx_t>( loc_col ) );
                coeffs_od.push_back( static_cast<scalar_t>( rvals[i] ) );
            }
        }
        nnz_d.push_back( nnzd );
        nnz_od.push_back( nnzod );
    }

    // Pad cols and cols_loc to ensure that all remote_DOFs are actually used
    nnz_pad = 0;
    if ( have_ghosts ) {
      std::set<gidx_t> colSet( cols_od.begin(), cols_od.end() );
      for ( auto rd : remote_DOFs ) {
	auto col = static_cast<gidx_t>( rd );
	auto cs = colSet.insert( col );
	if ( cs.second ) {
	  cols_od.push_back( col );
	  auto loc_col = comm_list->getLocalGhostID( rd );
	  cols_loc_od.push_back( static_cast<lidx_t>( loc_col ) );
	  coeffs_od.push_back( 0.0 );
	  nnz_pad++;
	}
      }
      nnz_od.back() += nnz_pad;
    }

    // Fill in row starts from nnz patterns
    const auto nRows = endRow - firstRow;
    rowstart_d.resize( nRows + 1 );
    rowstart_od.resize( nRows + 1 );
    std::exclusive_scan( nnz_d.begin(), nnz_d.end(), rowstart_d.begin(), 0 );
    std::exclusive_scan( nnz_od.begin(), nnz_od.end(), rowstart_od.begin(), 0 );
    rowstart_d[nRows] = rowstart_d[nRows - 1] + nnz_d[nRows - 1];
    rowstart_od[nRows] = rowstart_od[nRows - 1] + nnz_od[nRows - 1];
}
} // namespace AMP::LinearAlgebra

#endif
