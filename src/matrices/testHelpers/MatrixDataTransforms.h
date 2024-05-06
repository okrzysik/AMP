#ifndef included_MatrixDataTransforms_h
#define included_MatrixDataTransforms_h

#include "AMP/matrices/Matrix.h"

#include <type_traits>

namespace AMP::LinearAlgebra {

template<typename Policy>
void transformDofToCSR( std::shared_ptr<Matrix> matrix,
                        typename Policy::gidx_t &firstRow,
                        typename Policy::gidx_t &endRow,
                        std::vector<typename Policy::lidx_t> &nnz_d,
                        std::vector<typename Policy::gidx_t> &cols_d,
                        std::vector<typename Policy::lidx_t> &cols_loc_d,
                        std::vector<typename Policy::scalar_t> &coeffs_d,
                        std::vector<typename Policy::lidx_t> &nnz_od,
                        std::vector<typename Policy::gidx_t> &cols_od,
                        std::vector<typename Policy::lidx_t> &cols_loc_od,
                        std::vector<typename Policy::scalar_t> &coeffs_od )
{
    AMP_ASSERT( matrix );

    firstRow = static_cast<typename Policy::gidx_t>( matrix->beginRow() );
    endRow   = static_cast<typename Policy::gidx_t>( matrix->endRow() );

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
    for ( auto row = firstRow; row < endRow; ++row ) {

        std::vector<size_t> rcols;
        std::vector<typename Policy::scalar_t> rvals;

        matrix->getRowByGlobalID( row, rcols, rvals );

        // Loop over columns and insert into on and off diagonal blocks
        typename Policy::lidx_t nnzd = 0, nnzod = 0;
        for ( size_t i = 0; i < rcols.size(); ++i ) {
            const auto col = static_cast<typename Policy::gidx_t>( rcols[i] );
            if ( firstRow <= col && col < endRow ) {
                nnzd++;
                cols_d.push_back( col );
                cols_loc_d.push_back( static_cast<typename Policy::lidx_t>( col - firstRow ) );
                coeffs_d.push_back( static_cast<typename Policy::scalar_t>( rvals[i] ) );
            } else if ( have_ghosts ) {
                nnzod++;
                cols_od.push_back( col );
                auto loc_col = comm_list->getLocalGhostID( rcols[i] );
                cols_loc_od.push_back( static_cast<typename Policy::lidx_t>( loc_col ) );
                coeffs_od.push_back( static_cast<typename Policy::scalar_t>( rvals[i] ) );
            }
        }
        nnz_d.push_back( nnzd );
        nnz_od.push_back( nnzod );
    }
}
} // namespace AMP::LinearAlgebra

#endif
