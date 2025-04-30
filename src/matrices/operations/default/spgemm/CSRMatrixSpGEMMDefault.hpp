#include "AMP/matrices/operations/default/spgemm/CSRMatrixSpGEMMDefault.h"

#include "ProfilerApp.h"

#include <iostream>
#include <map>
#include <set>
#include <unordered_map>

namespace AMP::LinearAlgebra {

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::symbolicMultiply()
{
    PROFILE( "CSRMatrixSpGEMMDefault::symbolicMultiply" );

    // start communication to build BRemote before doing anything
    if ( A->hasOffDiag() ) {
        startBRemoteComm();
    }

    auto A_diag = A->getDiagMatrix();
    auto A_offd = A->getOffdMatrix();
    auto B_diag = B->getDiagMatrix();
    auto B_offd = B->getOffdMatrix();
    C_diag_diag = std::make_shared<DiagMatrixData>( nullptr,
                                                    C->getMemoryLocation(),
                                                    C->beginRow(),
                                                    C->endRow(),
                                                    C->beginCol(),
                                                    C->endCol(),
                                                    true );
    C_diag_offd = std::make_shared<DiagMatrixData>( nullptr,
                                                    C->getMemoryLocation(),
                                                    C->beginRow(),
                                                    C->endRow(),
                                                    C->beginCol(),
                                                    C->endCol(),
                                                    false );

    {
        PROFILE( "CSRMatrixSpGEMMDefault::symbolicMultiply (A_diag)" );
        multiply<DenseAccumulator, true>( A_diag, B_diag, C_diag_diag );
        multiply<DenseAccumulator, true>( A_diag, B_offd, C_diag_offd );
    }

    if ( A->hasOffDiag() ) {
        PROFILE( "CSRMatrixSpGEMMDefault::symbolicMultiply (A_offd)" );
        endBRemoteComm();
        C_offd_diag = std::make_shared<DiagMatrixData>( nullptr,
                                                        C->getMemoryLocation(),
                                                        C->beginRow(),
                                                        C->endRow(),
                                                        C->beginCol(),
                                                        C->endCol(),
                                                        true );
        C_offd_offd = std::make_shared<DiagMatrixData>( nullptr,
                                                        C->getMemoryLocation(),
                                                        C->beginRow(),
                                                        C->endRow(),
                                                        C->beginCol(),
                                                        C->endCol(),
                                                        false );

        multiply<DenseAccumulator, true>( A_offd, BR_diag, C_offd_diag );
        multiply<DenseAccumulator, true>( A_offd, BR_offd, C_offd_offd );
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::numericMultiply()
{
    PROFILE( "CSRMatrixSpGEMMDefault::numericMultiply" );

    // start communication to build BRemote before doing anything
    if ( A->hasOffDiag() && d_need_comms ) {
        startBRemoteComm();
    }

    auto A_diag = A->getDiagMatrix();
    auto A_offd = A->getOffdMatrix();
    auto B_diag = B->getDiagMatrix();
    auto B_offd = B->getOffdMatrix();
    auto C_diag = C->getDiagMatrix();
    auto C_offd = C->getOffdMatrix();

    {
        PROFILE( "CSRMatrixSpGEMMDefault::numericMultiply (A_diag)" );
        multiply<DenseAccumulator, false>( A_diag, B_diag, C_diag_diag );
        multiply<DenseAccumulator, false>( A_diag, B_offd, C_diag_offd );
    }

    if ( A->hasOffDiag() ) {
        PROFILE( "CSRMatrixSpGEMMDefault::numericMultiply (A_offd)" );
        if ( d_need_comms ) {
            endBRemoteComm();
        }
        multiply<DenseAccumulator, false>( A_offd, BR_diag, C_offd_diag );
        multiply<DenseAccumulator, false>( A_offd, BR_offd, C_offd_offd );
    }

    // merge the separate blocks together into cohesive output matrix
    // and deallocate blocks now that they are not needed
    C_diag->mergeMatrices( C_diag_diag, C_offd_diag );
    C_offd->mergeMatrices( C_diag_offd, C_offd_offd );
    C_diag_diag.reset();
    C_offd_diag.reset();
    C_diag_offd.reset();
    C_offd_offd.reset();

    C->globalToLocalColumns();
    C->resetDOFManagers();

    // set that comms need to be refreshed
    // assumes that user will only call multiply again if they have changed
    // the values in A and or B
    d_need_comms = true;
}

template<typename Policy, class Allocator, class DiagMatrixData>
template<class Accumulator, bool IsSymbolic>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::multiply(
    std::shared_ptr<DiagMatrixData> A_data,
    std::shared_ptr<DiagMatrixData> B_data,
    std::shared_ptr<DiagMatrixData> C_data )
{
    PROFILE( "CSRMatrixSpGEMMDefault::multiplyBlock" );

    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    AMP_DEBUG_ASSERT( A_data != nullptr );
    AMP_DEBUG_ASSERT( B_data != nullptr );
    AMP_DEBUG_ASSERT( C_data != nullptr );

    if ( A_data->isEmpty() || B_data->isEmpty() ) {
        return;
    }

    const auto nRows     = static_cast<lidx_t>( A->numLocalRows() );
    const bool B_is_diag = B_data->isDiag();

    // all fields from blocks involved
    lidx_t *A_rs = nullptr, *A_cols_loc = nullptr;
    gidx_t *A_cols     = nullptr;
    scalar_t *A_coeffs = nullptr;

    lidx_t *B_rs = nullptr, *B_cols_loc = nullptr;
    gidx_t *B_cols     = nullptr;
    scalar_t *B_coeffs = nullptr;

    lidx_t *C_rs = nullptr, *C_cols_loc = nullptr;
    gidx_t *C_cols     = nullptr;
    scalar_t *C_coeffs = nullptr;

    // Extract available fields
    std::tie( A_rs, A_cols, A_cols_loc, A_coeffs ) = A_data->getDataFields();
    std::tie( B_rs, B_cols, B_cols_loc, B_coeffs ) = B_data->getDataFields();
    std::tie( C_rs, C_cols, C_cols_loc, C_coeffs ) = C_data->getDataFields();

    AMP_DEBUG_ASSERT( A_cols_loc != nullptr );
    AMP_DEBUG_ASSERT( B_cols_loc != nullptr );

    // may or may not have access to B global column indices
    // set up conversion function from local indices
    auto B_colmap          = B_data->getColumnMap();
    auto B_colmap_size     = B_data->numUniqueColumns();
    const auto B_first_col = B_data->beginCol();

    auto B_to_global = [B_cols_loc, B_first_col, B_colmap, B_is_diag]( const lidx_t k ) -> gidx_t {
        return B_is_diag ? B_first_col + B_cols_loc[k] : B_colmap[B_cols_loc[k]];
    };

    // Create accumulator with appropriate capacity
    lidx_t acc_cap = B_is_diag ? B_data->numLocalColumns() : B_colmap_size;
    Accumulator acc( acc_cap );

    // Finally, after all the setup do the actual computation
    if constexpr ( IsSymbolic ) {
        // If this is a symbolic call just count NZ and write to
        // rs field in C
        for ( lidx_t row = 0; row < nRows; ++row ) {
            // get rows in B block from the A_diag column indices
            for ( lidx_t j = A_rs[row]; j < A_rs[row + 1]; ++j ) {
                const auto Acl = A_cols_loc[j];
                // then row of C is union of those B row nz patterns
                for ( lidx_t k = B_rs[Acl]; k < B_rs[Acl + 1]; ++k ) {
                    const auto bc = B_to_global( k );
                    acc.insert_or_append( B_cols_loc[k], bc );
                }
            }
            // write out row length and clear accumulator
            C_rs[row] += acc.num_inserted;
            acc.clear();
        }
        C_data->setNNZ( true );
    } else {
        // Otherwise, for numeric call write directly into C by
        // passing pointers into cols and coeffs fields as workspace
        // for the accumulator
        for ( lidx_t row = 0; row < nRows; ++row ) {
            const auto row_len = C_rs[row + 1] - C_rs[row];
            auto cols          = &C_cols[C_rs[row]];
            auto vals          = &C_coeffs[C_rs[row]];
            // get rows in B block from the A column indices
            for ( lidx_t j = A_rs[row]; j < A_rs[row + 1]; ++j ) {
                const auto Acl  = A_cols_loc[j];
                const auto Aval = A_coeffs[j];
                // then row of C is union of those B row nz patterns
                for ( lidx_t k = B_rs[Acl]; k < B_rs[Acl + 1]; ++k ) {
                    const auto bc = B_to_global( k );
                    acc.insert_or_append(
                        B_cols_loc[k], bc, Aval * B_coeffs[k], cols, vals, row_len );
                }
            }
            acc.clear();
        }
        C_data->sortColumns();
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::setupBRemoteComm()
{
    /*
     * Setting up the comms is somewhat involved. A high level overview
     * of the steps involved is:
     * 1. Collect comm list info and needed remote rows
     * 2. Trim down lists from steps 3 and 4 to ranks that are actually needed
     * 3. Record which specific rows are needed from each process
     * 4. Send row ids from 6 to owners of those rows

     * NOTES:
     *  Step 4 uses non-blocking recvs and blocking sends.
     */

    PROFILE( "CSRMatrixSpGEMMDefault::setupBRemoteComm" );

    using lidx_t = typename Policy::lidx_t;

    auto comm_size = comm.getSize();

    // 1. Query comm list info and get offd colmap
    auto comm_list            = A->getRightCommList();
    auto rows_per_rank_recv   = comm_list->getReceiveSizes();
    auto rows_per_rank_send   = comm_list->getSendSizes();
    auto B_last_rows          = comm_list->getPartition();
    const auto A_col_map_size = A->getOffdMatrix()->numUniqueColumns();
    const auto A_col_map      = A->getOffdMatrix()->getColumnMap();

    // 2. the above rows per rank lists generally include lots of zeros
    // trim down to the ranks that actually need to communicate
    int total_send = 0, total_recv = 0;
    for ( int r = 0; r < comm_size; ++r ) {
        const auto nsend = rows_per_rank_send[r];
        if ( nsend > 0 ) {
            d_dest_info.insert( std::make_pair( r, SpGEMMCommInfo( nsend ) ) );
        }
        const auto nrecv = rows_per_rank_recv[r];
        if ( nrecv > 0 ) {
            d_src_info.insert( std::make_pair( r, SpGEMMCommInfo( nrecv ) ) );
        }
        total_send += nsend;
        total_recv += nrecv;
    }

    // 3. Scan over column map now writing into the trimmed down src list
    for ( lidx_t n = 0; n < A_col_map_size; ++n ) {
        const auto col = static_cast<std::size_t>( A_col_map[n] );
        int owner      = -1;
        if ( col < B_last_rows[0] ) {
            owner = 0;
        } else {
            for ( int r = 1; r < comm_size; ++r ) {
                auto rs = B_last_rows[r - 1], re = B_last_rows[r];
                if ( rs <= col && col < re ) {
                    owner = r;
                    break;
                }
            }
        }
        d_src_info[owner].rowids.push_back( col );
    }

    // 4. send rowids to their owners
    // start by posting the irecvs
    const int TAG = 7800;
    std::vector<AMP_MPI::Request> irecvs;
    for ( auto it = d_dest_info.begin(); it != d_dest_info.end(); ++it ) {
        it->second.rowids.resize( it->second.numrow );
        irecvs.push_back(
            comm.Irecv( it->second.rowids.data(), it->second.numrow, it->first, TAG ) );
    }
    // now send all the rows we want from other ranks
    for ( auto it = d_src_info.begin(); it != d_src_info.end(); ++it ) {
        comm.send( it->second.rowids.data(), it->second.numrow, it->first, TAG );
    }
    // wait for receives to finish
    comm.waitAll( static_cast<int>( irecvs.size() ), irecvs.data() );
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::startBRemoteComm()
{
    // check if the communicator information is available and create if needed
    if ( d_dest_info.empty() ) {
        setupBRemoteComm();
    }

    // subset matrices by rows that other ranks need and send them out
    for ( auto it = d_dest_info.begin(); it != d_dest_info.end(); ++it ) {
        auto block = B->subsetRows( it->second.rowids );
        d_send_matrices.insert( { it->first, block } );
    }
    d_csr_comm.sendMatrices( d_send_matrices );
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::endBRemoteComm()
{
    using lidx_t = typename Policy::lidx_t;

    PROFILE( "CSRMatrixSpGEMMDefault::endBRemoteComm" );

    d_recv_matrices = d_csr_comm.recvMatrices( 0, 0, 0, B->numGlobalColumns() );
    // BRemotes do not need any particular parameters object internally
    BR_diag = CSRLocalMatrixData<Policy, Allocator>::ConcatVertical(
        nullptr, d_recv_matrices, B->beginCol(), B->endCol(), true );
    BR_offd = CSRLocalMatrixData<Policy, Allocator>::ConcatVertical(
        nullptr, d_recv_matrices, B->beginCol(), B->endCol(), false );
    const auto A_col_map_size = A->getOffdMatrix()->numUniqueColumns();
    if ( A_col_map_size != static_cast<lidx_t>( BR_diag->endRow() ) ) {
        int num_reqd = 0;
        for ( auto it = d_src_info.begin(); it != d_src_info.end(); ++it ) {
            num_reqd += it->second.numrow;
        }
        std::cout << "Rank " << comm.getRank() << " expected last row " << A_col_map_size << " got "
                  << BR_diag->endRow() << " requested " << num_reqd << std::endl;

        AMP_ERROR( "BRemote has wrong ending row" );
    }

    // comms are done and BR_{diag,offd} filled, deallocate send/recv blocks
    d_send_matrices.clear();
    d_recv_matrices.clear();

    // set flag that recv'd matrices are valid
    d_need_comms = false;
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::DenseAccumulator::
    insert_or_append( typename Policy::lidx_t loc, typename Policy::gidx_t gbl )
{
    using lidx_t = typename Policy::lidx_t;

    const auto k = flags[loc];
    if ( k == -1 ) {
        flags[loc] = num_inserted;
        if ( num_inserted == static_cast<lidx_t>( flag_inv.size() ) ) {
            flag_inv.push_back( loc );
            cols.push_back( gbl );
        } else {
            flag_inv[num_inserted] = loc;
            cols[num_inserted]     = gbl;
        }
        ++num_inserted;
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::DenseAccumulator::
    insert_or_append( typename Policy::lidx_t loc,
                      typename Policy::gidx_t gbl,
                      typename Policy::scalar_t val,
                      typename Policy::gidx_t *col_space,
                      typename Policy::scalar_t *val_space,
                      [[maybe_unused]] typename Policy::lidx_t max_pos )
{
    using lidx_t = typename Policy::lidx_t;

    const auto k = flags[loc];
    if ( k == -1 ) {
        AMP_DEBUG_ASSERT( num_inserted < max_pos );
        flags[loc] = num_inserted;
        if ( num_inserted == static_cast<lidx_t>( flag_inv.size() ) ) {
            flag_inv.push_back( loc );
        } else {
            flag_inv[num_inserted] = loc;
        }
        col_space[num_inserted] = gbl;
        val_space[num_inserted] = val;
        ++num_inserted;
    } else {
        val_space[k] += val;
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::DenseAccumulator::clear()
{
    using lidx_t = typename Policy::lidx_t;

    for ( lidx_t n = 0; n < num_inserted; ++n ) {
        flags[flag_inv[n]] = -1;
    }
    num_inserted = 0;
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::SparseAccumulator::
    insert_or_append( typename Policy::lidx_t loc, typename Policy::gidx_t gbl )
{
    kv[gbl]      = loc;
    num_inserted = static_cast<lidx_t>( kv.size() );
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::SparseAccumulator::
    insert_or_append( [[maybe_unused]] typename Policy::lidx_t loc,
                      typename Policy::gidx_t gbl,
                      typename Policy::scalar_t val,
                      typename Policy::gidx_t *col_space,
                      typename Policy::scalar_t *val_space,
                      [[maybe_unused]] typename Policy::lidx_t max_pos )
{
    auto it = kv.find( gbl );
    if ( it != kv.end() ) {
        val_space[it->second] += val;
    } else {
        AMP_DEBUG_ASSERT( num_inserted < max_pos );
        kv.insert( std::make_pair( gbl, num_inserted ) );
        col_space[num_inserted] = gbl;
        val_space[num_inserted] = val;
    }
    num_inserted = static_cast<lidx_t>( kv.size() );
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::SparseAccumulator::clear()
{
    num_inserted = 0;
    kv.clear();
}

} // namespace AMP::LinearAlgebra
