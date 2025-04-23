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

    using lidx_t = typename Policy::lidx_t;

    // start communication to build BRemote before doing anything
    if ( A->hasOffDiag() ) {
        startBRemoteComm();
        endBRemoteComm();
    }

    const auto nRows = static_cast<lidx_t>( A->numLocalRows() );

    auto B_diag = B->getDiagMatrix();
    auto B_offd = B->getOffdMatrix();

    std::vector<lidx_t> nnz_diag( nRows, 0 ), nnz_offd( nRows, 0 );

    multiplyLocal<true>( B_diag, nullptr, nnz_diag.data() );
    if ( B->hasOffDiag() ) {
        multiplyLocal<true>( B_offd, nullptr, nnz_offd.data() );
    }

    // Give C the nnz counts so that it can allocate space internally
    C->setNNZ( nnz_diag, nnz_offd );
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::numericMultiply()
{
    PROFILE( "CSRMatrixSpGEMMDefault::numericMultiply" );

    // start communication to build BRemote before doing anything
    if ( A->hasOffDiag() && d_need_comms ) {
        startBRemoteComm();
        endBRemoteComm();
    }

    auto B_diag = B->getDiagMatrix();
    auto B_offd = B->getOffdMatrix();
    auto C_diag = C->getDiagMatrix();
    auto C_offd = C->getOffdMatrix();

    // Process diagonal block of A acting on whole local part of B
    multiplyLocal<false>( B_diag, C_diag, nullptr );
    if ( B->hasOffDiag() ) {
        multiplyLocal<false>( B_offd, C_offd, nullptr );
    }

    C->globalToLocalColumns();
    C->resetDOFManagers();

    // set that comms need to be refreshed
    // assumes that user will only call multiply again if they have changed
    // the values in A and or B
    d_need_comms = true;
}

template<typename Policy, class Allocator, class DiagMatrixData>
template<bool SYMBOLIC>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::multiplyLocal(
    std::shared_ptr<DiagMatrixData> B_data,
    std::shared_ptr<DiagMatrixData> C_data,
    typename Policy::lidx_t *nnz )
{
    PROFILE( "CSRMatrixSpGEMMDefault::multiplyLocal" );

    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    // need both blocks from A
    auto A_diag            = A->getDiagMatrix();
    auto A_offd            = A->getOffdMatrix();
    const auto A_have_offd = A->hasOffDiag();
    const auto nRows       = static_cast<lidx_t>( A->numLocalRows() );

    // all fields from blocks involved
    lidx_t *A_rs_d = nullptr, *A_cols_loc_d = nullptr;
    gidx_t *A_cols_d     = nullptr;
    scalar_t *A_coeffs_d = nullptr;

    lidx_t *A_rs_od = nullptr, *A_cols_loc_od = nullptr;
    gidx_t *A_cols_od     = nullptr;
    scalar_t *A_coeffs_od = nullptr;

    lidx_t *B_rs = nullptr, *B_cols_loc = nullptr;
    gidx_t *B_cols     = nullptr;
    scalar_t *B_coeffs = nullptr;

    lidx_t *BR_rs = nullptr, *BR_cols_loc = nullptr;
    gidx_t *BR_cols     = nullptr;
    scalar_t *BR_coeffs = nullptr;

    lidx_t *C_rs = nullptr, *C_cols_loc = nullptr;
    gidx_t *C_cols     = nullptr;
    scalar_t *C_coeffs = nullptr;

    // Extract available fields
    std::tie( A_rs_d, A_cols_d, A_cols_loc_d, A_coeffs_d ) = A_diag->getDataFields();
    std::tie( B_rs, B_cols, B_cols_loc, B_coeffs )         = B_data->getDataFields();
    if ( A_have_offd ) {
        AMP_DEBUG_ASSERT( BRemote != nullptr );
        std::tie( A_rs_od, A_cols_od, A_cols_loc_od, A_coeffs_od ) = A_offd->getDataFields();
        std::tie( BR_rs, BR_cols, BR_cols_loc, BR_coeffs )         = BRemote->getDataFields();
    }
    if constexpr ( !SYMBOLIC ) {
        AMP_DEBUG_ASSERT( C_data != nullptr );
        std::tie( C_rs, C_cols, C_cols_loc, C_coeffs ) = C_data->getDataFields();
    }

    // The output is for a block of C matching the type of the given B block
    const bool B_is_diag = B_data->isDiag();

    // Column range of diagonal C block is row range of B
    const auto col_start = B_data->beginRow();
    const auto col_end   = B_data->endRow();

    // test if a given global ID lands in/out of column range
    // depending on diag vs offd status
    auto idx_test = [col_start, col_end, B_is_diag]( const gidx_t col ) -> bool {
        return B_is_diag ? ( col_start <= col && col < col_end ) :
                           ( col < col_start || col_end <= col );
    };

    // may or may not have access to B global column indices
    // set up conversion function from local indices
    auto B_colmap          = B_data->getColumnMap();
    auto B_colmap_size     = B_data->numUniqueColumns();
    const auto B_first_col = B_data->beginCol();
    const bool have_B_cols = ( B_cols != nullptr );

    auto B_to_global = [B_cols, B_cols_loc, B_first_col, B_colmap, B_is_diag, have_B_cols](
                           const lidx_t k ) -> gidx_t {
        return have_B_cols ? B_cols[k] :
                             ( B_is_diag ? B_first_col + B_cols_loc[k] : B_colmap[B_cols_loc[k]] );
    };

    // In the opposite direction, BRemote does not have local indices
    // create similar conversion. This needs a unified column map over
    // both B_offd and BRemote when constructing C_offd.
    std::unordered_map<gidx_t, lidx_t> joined_b_colmap;
    if ( A_have_offd ) {
        for ( lidx_t loc = 0; loc < B_colmap_size; ++loc ) {
            joined_b_colmap[B_colmap[loc]] = loc;
        }
        for ( lidx_t n = 0; n < BRemote->numberOfNonZeros(); ++n ) {
            joined_b_colmap.insert(
                std::make_pair( BR_cols[n], static_cast<lidx_t>( joined_b_colmap.size() ) ) );
        }
    }
    auto BR_to_local =
        [B_first_col, &joined_b_colmap, B_is_diag, B_colmap_size]( const gidx_t k ) -> lidx_t {
        if ( B_is_diag ) {
            return static_cast<lidx_t>( k - B_first_col );
        } else {
            // this branch is expensive, though it is only called
            // for positions coming from A_offd * BRemote_offd
            // which is hopefully rare relatively speaking...
            auto it = joined_b_colmap.find( k );
            if ( it != joined_b_colmap.end() ) {
                return it->second;
            }
        }
        AMP_ERROR( "CSRMatrixSpGEMMDefault::symbolicMultiplyLocal BR_to_local failed" );
        return -1;
    };

    // create a dense accumulator based on the shape of B
    // rows in C are unions of rows of B, so can't be larger than
    // B span
    const lidx_t acc_capacity = B_is_diag ? B_data->numLocalColumns() : joined_b_colmap.size();
    DenseAccumulator acc( acc_capacity );

    // for each row in A
    if constexpr ( SYMBOLIC ) {
        for ( lidx_t row = 0; row < nRows; ++row ) {
            // get rows in B block from the A_diag column indices
            for ( lidx_t j = A_rs_d[row]; j < A_rs_d[row + 1]; ++j ) {
                const auto Acl = A_cols_loc_d[j];
                // then row of C is union of those B row nz patterns
                for ( lidx_t k = B_rs[Acl]; k < B_rs[Acl + 1]; ++k ) {
                    const auto bc = B_to_global( k );
                    if ( idx_test( bc ) ) {
                        acc.insert_or_append( B_cols_loc[k], bc );
                    }
                }
            }

            // do same for A_offd acting on BRemote if needed
            if ( A_have_offd ) {
                for ( lidx_t j = A_rs_od[row]; j < A_rs_od[row + 1]; ++j ) {
                    const auto Acl = A_cols_loc_od[j];
                    // then row of C is union of those B row nz patterns
                    for ( lidx_t k = BR_rs[Acl]; k < BR_rs[Acl + 1]; ++k ) {
                        const auto bc = BR_cols[k];
                        if ( idx_test( bc ) ) {
                            const auto loc = BR_to_local( bc );
                            acc.insert_or_append( loc, bc );
                        }
                    }
                }
            }

            // write out row length and clear accumulator
            nnz[row] = acc.num_inserted;
            acc.clear();
        }
    } else {
        // for each row in A block
        for ( lidx_t row = 0; row < nRows; ++row ) {
            auto cols = &C_cols[C_rs[row]];
            auto vals = &C_coeffs[C_rs[row]];
            // get rows in B block from the A column indices
            for ( lidx_t j = A_rs_d[row]; j < A_rs_d[row + 1]; ++j ) {
                const auto Acl  = A_cols_loc_d[j];
                const auto Aval = A_coeffs_d[j];
                // then row of C is union of those B row nz patterns
                for ( lidx_t k = B_rs[Acl]; k < B_rs[Acl + 1]; ++k ) {
                    const auto bc = B_to_global( k );
                    if ( idx_test( bc ) ) {
                        acc.insert_or_append( B_cols_loc[k], bc, Aval * B_coeffs[k], cols, vals );
                    }
                }
            }

            // do same for A_offd acting on BRemote if needed
            if ( A_have_offd ) {
                for ( lidx_t j = A_rs_od[row]; j < A_rs_od[row + 1]; ++j ) {
                    const auto Acl  = A_cols_loc_od[j];
                    const auto Aval = A_coeffs_od[j];
                    // then row of C is union of those B row nz patterns
                    for ( lidx_t k = BR_rs[Acl]; k < BR_rs[Acl + 1]; ++k ) {
                        const auto bc = BR_cols[k];
                        if ( idx_test( bc ) ) {
                            const auto loc = BR_to_local( bc );
                            acc.insert_or_append( loc, bc, Aval * BR_coeffs[k], cols, vals );
                        }
                    }
                }
            }
            acc.clear();
        }
    }
}

// template<typename Policy, class Allocator, class DiagMatrixData>
// void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::symbolicMultiplyLocal(
//     std::shared_ptr<DiagMatrixData> B_data, std::vector<typename Policy::lidx_t> &nnz )
// {
//     PROFILE( "CSRMatrixSpGEMMDefault::symbolicMultiplyLocal" );

//     using lidx_t   = typename Policy::lidx_t;
//     using gidx_t   = typename Policy::gidx_t;
//     using scalar_t = typename Policy::scalar_t;

//     // need both blocks from A
//     auto A_diag            = A->getDiagMatrix();
//     auto A_offd            = A->getOffdMatrix();
//     const auto A_have_offd = A->hasOffDiag();
//     const auto nRows       = static_cast<lidx_t>( A->numLocalRows() );

//     // all fields from blocks involved
//     lidx_t *A_rs_d = nullptr, *A_cols_loc_d = nullptr;
//     gidx_t *A_cols_d     = nullptr;
//     scalar_t *A_coeffs_d = nullptr;

//     lidx_t *A_rs_od = nullptr, *A_cols_loc_od = nullptr;
//     gidx_t *A_cols_od     = nullptr;
//     scalar_t *A_coeffs_od = nullptr;

//     lidx_t *B_rs = nullptr, *B_cols_loc = nullptr;
//     gidx_t *B_cols     = nullptr;
//     scalar_t *B_coeffs = nullptr;

//     lidx_t *BR_rs = nullptr, *BR_cols_loc = nullptr;
//     gidx_t *BR_cols     = nullptr;
//     scalar_t *BR_coeffs = nullptr;

//     // Extract available fields
//     std::tie( A_rs_d, A_cols_d, A_cols_loc_d, A_coeffs_d ) = A_diag->getDataFields();
//     std::tie( B_rs, B_cols, B_cols_loc, B_coeffs )         = B_data->getDataFields();
//     if ( A_have_offd ) {
//         AMP_DEBUG_ASSERT( BRemote != nullptr );
//         std::tie( A_rs_od, A_cols_od, A_cols_loc_od, A_coeffs_od ) = A_offd->getDataFields();
//         std::tie( BR_rs, BR_cols, BR_cols_loc, BR_coeffs )         = BRemote->getDataFields();
//     }

//     // The output is for a block of C matching the type of the given B block
//     const bool B_is_diag = B_data->isDiag();

//     // Column range of diagonal C block is row range of B
//     const auto col_start = B_data->beginRow();
//     const auto col_end   = B_data->endRow();

//     // test if a given global ID lands in/out of column range
//     // depending on diag vs offd status
//     auto idx_test = [col_start, col_end, B_is_diag]( const gidx_t col ) -> bool {
//         return B_is_diag ? ( col_start <= col && col < col_end ) :
//                            ( col < col_start || col_end <= col );
//     };

//     // create a dense accumulator based on the shape of B
//     // rows in C are unions of rows of B, so can't be larger than
//     // B span
//     DenseAccumulator acc( B_is_diag ? B_data->numLocalColumns() : B_data->numUniqueColumns() );

//     // may or may not have access to B global column indices
//     // set up conversion function from local indices
//     auto B_colmap          = B_data->getColumnMap();
//     auto B_colmap_size     = B_data->numUniqueColumns();
//     const auto B_first_col = B_data->beginCol();
//     const bool have_B_cols = ( B_cols != nullptr );

//     auto B_to_global = [B_cols, B_cols_loc, B_first_col, B_colmap, B_is_diag, have_B_cols](
//                            const lidx_t k ) -> gidx_t {
//         return have_B_cols ? B_cols[k] :
//                              ( B_is_diag ? B_first_col + B_cols_loc[k] : B_colmap[B_cols_loc[k]]
//                              );
//     };

//     // In the opposite direction, BRemote does not have local indices
//     // create similar conversion
//     // ******** This does not work... Need single colmap for both B and BR ********
//     auto BR_to_local =
//         [B_first_col, B_colmap, B_colmap_size, B_is_diag]( const gidx_t k ) -> lidx_t {
//         if ( B_is_diag ) {
//             return static_cast<lidx_t>( k - B_first_col );
//         } else {
//             // this branch is expensive, though it is only called
//             // for positions coming from A_offd * BRemote_offd
//             // which is hopefully rare relatively speaking...
//             for ( lidx_t loc = 0; loc < B_colmap_size; ++loc ) {
//                 if ( B_colmap[loc] == k ) {
//                     return loc;
//                 }
//             }
//         }
//         AMP_ERROR( "CSRMatrixSpGEMMDefault::symbolicMultiplyLocal BR_to_local failed" );
//         return -1;
//     };

//     // for each row in A
//     for ( lidx_t row = 0; row < nRows; ++row ) {
//         // get rows in B block from the A_diag column indices
//         for ( lidx_t j = A_rs_d[row]; j < A_rs_d[row + 1]; ++j ) {
//             auto Acl = A_cols_loc_d[j];
//             // then row of C is union of those B row nz patterns
//             for ( lidx_t k = B_rs[Acl]; k < B_rs[Acl + 1]; ++k ) {
//                 const auto bc = B_to_global( k );
//                 if ( idx_test( bc ) ) {
//                     acc.insert_or_append( B_cols_loc[k], bc );
//                 }
//             }
//         }

//         // do same for A_offd acting on BRemote if needed
//         if ( A_have_offd ) {
//             for ( lidx_t j = A_rs_od[row]; j < A_rs_od[row + 1]; ++j ) {
//                 auto Acl = A_cols_loc_od[j];
//                 // then row of C is union of those B row nz patterns
//                 for ( lidx_t k = BR_rs[Acl]; k < BR_rs[Acl + 1]; ++k ) {
//                     const auto bc = BR_cols[k];
//                     if ( idx_test( bc ) ) {
//                         const auto loc = BR_to_local( bc );
//                         acc.insert_or_append( loc, bc );
//                     }
//                 }
//             }
//         }

//         // write out row length and clear accumulator
//         nnz[row] = acc.num_inserted;
//         acc.clear();
//     }
// }

// template<typename Policy, class Allocator, class DiagMatrixData>
// void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::numericMultiplyLocal(
//     std::shared_ptr<DiagMatrixData> B_data, std::shared_ptr<DiagMatrixData> C_data )
// {
//     PROFILE( "CSRMatrixSpGEMMDefault::numericMultiplyLocal" );

//     using lidx_t   = typename Policy::lidx_t;
//     using gidx_t   = typename Policy::gidx_t;
//     using scalar_t = typename Policy::scalar_t;

//     // need both blocks from A
//     auto A_diag            = A->getDiagMatrix();
//     auto A_offd            = A->getOffdMatrix();
//     const auto A_have_offd = A->hasOffDiag();
//     const auto nRows       = static_cast<lidx_t>( A->numLocalRows() );

//     // all fields from blocks involved
//     lidx_t *A_rs_d = nullptr, *A_cols_loc_d = nullptr;
//     gidx_t *A_cols_d     = nullptr;
//     scalar_t *A_coeffs_d = nullptr;

//     lidx_t *A_rs_od = nullptr, *A_cols_loc_od = nullptr;
//     gidx_t *A_cols_od     = nullptr;
//     scalar_t *A_coeffs_od = nullptr;

//     lidx_t *B_rs = nullptr, *B_cols_loc = nullptr;
//     gidx_t *B_cols     = nullptr;
//     scalar_t *B_coeffs = nullptr;

//     lidx_t *BR_rs = nullptr, *BR_cols_loc = nullptr;
//     gidx_t *BR_cols     = nullptr;
//     scalar_t *BR_coeffs = nullptr;

//     lidx_t *C_rs = nullptr, *C_cols_loc = nullptr;
//     gidx_t *C_cols     = nullptr;
//     scalar_t *C_coeffs = nullptr;

//     // Extract available fields
//     std::tie( A_rs_d, A_cols_d, A_cols_loc_d, A_coeffs_d ) = A_diag->getDataFields();
//     std::tie( B_rs, B_cols, B_cols_loc, B_coeffs )         = B_data->getDataFields();
//     std::tie( C_rs, C_cols, C_cols_loc, C_coeffs )         = C_data->getDataFields();
//     if ( A_have_offd ) {
//         AMP_DEBUG_ASSERT( BRemote != nullptr );
//         std::tie( A_rs_od, A_cols_od, A_cols_loc_od, A_coeffs_od ) = A_offd->getDataFields();
//         std::tie( BR_rs, BR_cols, BR_cols_loc, BR_coeffs )         = BRemote->getDataFields();
//     }

//     // The output is for a block of C matching the type of the given B block
//     const bool B_is_diag = B_data->isDiag();

//     // Column range of diagonal C block is row range of B
//     const auto col_start = B_data->beginRow();
//     const auto col_end   = B_data->endRow();

//     // test if a given global ID lands in/out of column range
//     // depending on diag vs offd status
//     auto idx_test = [col_start, col_end, B_is_diag]( const gidx_t col ) -> bool {
//         return B_is_diag ? ( col_start <= col && col < col_end ) :
//                            ( col < col_start || col_end <= col );
//     };

//     // create a dense accumulator based on the shape of B
//     // rows in C are unions of rows of B, so can't be larger than
//     // B span
//     DenseAccumulator acc( B_is_diag ? B_data->numLocalColumns() : B_data->numUniqueColumns() );

//     // may or may not have access to B global column indices
//     // set up conversion function from local indices
//     auto B_colmap          = B_data->getColumnMap();
//     auto B_colmap_size     = B_data->numUniqueColumns();
//     const auto B_first_col = B_data->beginCol();
//     const bool have_B_cols = ( B_cols != nullptr );

//     auto B_to_global = [B_cols, B_cols_loc, B_first_col, B_colmap, B_is_diag, have_B_cols](
//                            const lidx_t k ) -> gidx_t {
//         return have_B_cols ? B_cols[k] :
//                              ( B_is_diag ? B_first_col + B_cols_loc[k] : B_colmap[B_cols_loc[k]]
//                              );
//     };

//     // In the opposite direction, BRemote does not have local indices
//     // create similar conversion
//     auto BR_to_local =
//         [B_first_col, B_colmap, B_colmap_size, B_is_diag]( const gidx_t k ) -> lidx_t {
//         if ( B_is_diag ) {
//             return static_cast<lidx_t>( k - B_first_col );
//         } else {
//             // this branch is expensive, though it is only called
//             // for positions coming from A_offd * BRemote_offd
//             // which is hopefully rare relatively speaking...
//             for ( lidx_t loc = 0; loc < B_colmap_size; ++loc ) {
//                 if ( B_colmap[loc] == k ) {
//                     return loc;
//                 }
//             }
//         }
//         AMP_ERROR( "CSRMatrixSpGEMMDefault::symbolicMultiplyLocal BR_to_local failed" );
//         return -1;
//     };

//     // for each row in A block
//     for ( lidx_t row = 0; row < nRows; ++row ) {
//         auto cols = &C_cols[C_rs[row]];
//         auto vals = &C_coeffs[C_rs[row]];
//         // get rows in B block from the A column indices
//         for ( lidx_t j = A_rs_d[row]; j < A_rs_d[row + 1]; ++j ) {
//             const auto Acl  = A_cols_loc_d[j];
//             const auto Aval = A_coeffs_d[j];
//             // then row of C is union of those B row nz patterns
//             for ( lidx_t k = B_rs[Acl]; k < B_rs[Acl + 1]; ++k ) {
//                 const auto bc = B_to_global( k );
//                 if ( idx_test( bc ) ) {
//                     acc.insert_or_append( B_cols_loc[k], bc, Aval * B_coeffs[k], cols, vals );
//                 }
//             }
//         }

//         // do same for A_offd acting on BRemote if needed
//         if ( A_have_offd ) {
//             for ( lidx_t j = A_rs_od[row]; j < A_rs_od[row + 1]; ++j ) {
//                 auto Acl        = A_cols_loc_od[j];
//                 const auto Aval = A_coeffs_od[j];
//                 // then row of C is union of those B row nz patterns
//                 for ( lidx_t k = BR_rs[Acl]; k < BR_rs[Acl + 1]; ++k ) {
//                     const auto bc = BR_cols[k];
//                     if ( idx_test( bc ) ) {
//                         const auto loc = BR_to_local( bc );
//                         acc.insert_or_append( loc, bc, Aval * B_coeffs[k], cols, vals );
//                     }
//                 }
//             }
//         }
//         acc.clear();
//     }
// }

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
    // BRemote does not need any particular parameters object internally
    BRemote = CSRLocalMatrixData<Policy, Allocator>::ConcatVertical( nullptr, d_recv_matrices );
    const auto A_col_map_size = A->getOffdMatrix()->numUniqueColumns();
    if ( A_col_map_size != static_cast<lidx_t>( BRemote->endRow() ) ) {
        int num_reqd = 0;
        for ( auto it = d_src_info.begin(); it != d_src_info.end(); ++it ) {
            num_reqd += it->second.numrow;
        }
        std::cout << "Rank " << comm.getRank() << " expected last row " << A_col_map_size << " got "
                  << BRemote->endRow() << " requested " << num_reqd << std::endl;

        AMP_ERROR( "BRemote has wrong ending row" );
    }

    // set flag that recv'd matrices are valid
    d_need_comms = false;
}
} // namespace AMP::LinearAlgebra
