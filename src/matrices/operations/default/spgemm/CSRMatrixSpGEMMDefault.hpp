#include "AMP/matrices/operations/default/spgemm/CSRMatrixSpGEMMDefault.h"

#include "ProfilerApp.h"

#ifndef CSRSPGEMM_REPORT_SPACC_STATS
    #define CSRSPGEMM_REPORT_SPACC_STATS 0
#endif

namespace AMP::LinearAlgebra {

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::symbolicMultiply()
{
    if ( d_overlap_comms ) {
        symbolicMultiply_Overlapped();
    } else {
        symbolicMultiply_NonOverlapped();
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::numericMultiply()
{
    if ( d_overlap_comms ) {
        numericMultiply_Overlapped();
    } else {
        numericMultiply_NonOverlapped();
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::
    symbolicMultiply_NonOverlapped()
{
    PROFILE( "symbolicMultiply_NonOverlapped" );

    // non-overlapped, so do full comms first
    if ( A->hasOffDiag() ) {
        startBRemoteComm();
        endBRemoteComm();
    }

    if ( !A->hasOffDiag() ) {
        multiply<DenseAccumulator, true>( A_diag, B_diag, C_diag );
        multiply<SparseAccumulator, true>( A_diag, B_offd, C_offd );
    } else {
        if ( BR_diag.get() ) {
            multiplyFused<DenseAccumulator, true, true>( B_diag, BR_diag, C_diag );
        } else {
            multiply<DenseAccumulator, true>( A_diag, B_diag, C_diag );
        }
        if ( BR_offd.get() ) {
            multiplyFused<SparseAccumulator, false, true>( B_offd, BR_offd, C_offd );
        } else {
            multiply<SparseAccumulator, true>( A_diag, B_offd, C_offd );
        }
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::symbolicMultiply_Overlapped()
{
    PROFILE( "symbolicMultiply_Overlapped" );

    // start communication to build BRemote before doing anything
    if ( A->hasOffDiag() ) {
        startBRemoteComm();
    }

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
        PROFILE( "symbolicMultiply_Overlapped (local)" );
        multiply<DenseAccumulator, true>( A_diag, B_diag, C_diag_diag );
        multiply<SparseAccumulator, true>( A_diag, B_offd, C_diag_offd );
    }

    if ( A->hasOffDiag() ) {
        PROFILE( "symbolicMultiply_Overlapped (remote)" );
        endBRemoteComm();
        if ( BR_diag.get() != nullptr ) {
            C_offd_diag = std::make_shared<DiagMatrixData>( nullptr,
                                                            C->getMemoryLocation(),
                                                            C->beginRow(),
                                                            C->endRow(),
                                                            C->beginCol(),
                                                            C->endCol(),
                                                            true );
            multiply<DenseAccumulator, true>( A_offd, BR_diag, C_offd_diag );
        }
        if ( BR_offd.get() != nullptr ) {
            C_offd_offd = std::make_shared<DiagMatrixData>( nullptr,
                                                            C->getMemoryLocation(),
                                                            C->beginRow(),
                                                            C->endRow(),
                                                            C->beginCol(),
                                                            C->endCol(),
                                                            false );
            multiply<SparseAccumulator, true>( A_offd, BR_offd, C_offd_offd );
        }
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::
    numericMultiply_NonOverlapped()
{
    PROFILE( "numericMultiply_NonOverlapped" );

    // non-overlapped comms so do full set first if needed
    if ( A->hasOffDiag() && d_need_comms ) {
        startBRemoteComm();
        endBRemoteComm();
    }

    if ( !A->hasOffDiag() ) {
        multiply<DenseAccumulator, false>( A_diag, B_diag, C_diag );
        multiply<SparseAccumulator, false>( A_diag, B_offd, C_offd );
    } else {
        if ( BR_diag.get() ) {
            multiplyFused<DenseAccumulator, true, false>( B_diag, BR_diag, C_diag );
        } else {
            multiply<DenseAccumulator, false>( A_diag, B_diag, C_diag );
        }
        if ( BR_offd.get() ) {
            multiplyFused<SparseAccumulator, false, false>( B_offd, BR_offd, C_offd );
        } else {
            multiply<SparseAccumulator, false>( A_diag, B_offd, C_offd );
        }
    }

    C->globalToLocalColumns();
    C->resetDOFManagers();

    // set that comms need to be refreshed
    // assumes that user will only call multiply again if they have changed
    // the values in A and or B
    d_need_comms = true;
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::numericMultiply_Overlapped()
{
    PROFILE( "numericMultiply_Overlapped" );

    // start communication to build BRemote before doing anything
    if ( A->hasOffDiag() && d_need_comms ) {
        startBRemoteComm();
    }

    {
        PROFILE( "numericMultiply_Overlapped (local)" );
        multiply<DenseAccumulator, false>( A_diag, B_diag, C_diag_diag );
        multiply<SparseAccumulator, false>( A_diag, B_offd, C_diag_offd );
    }

    if ( A->hasOffDiag() ) {
        PROFILE( "numericMultiply_Overlapped (remote)" );
        if ( d_need_comms ) {
            endBRemoteComm();
        }
        if ( BR_diag.get() != nullptr ) {
            multiply<DenseAccumulator, false>( A_offd, BR_diag, C_offd_diag );
        }
        if ( BR_offd.get() != nullptr ) {
            multiply<SparseAccumulator, false>( A_offd, BR_offd, C_offd_offd );
        }
    }

    {
        PROFILE( "numericMultiply_Overlapped (merge)" );
        mergeDiag<DenseAccumulator>();
        mergeOffd<SparseAccumulator>();
    }

    C->globalToLocalColumns();
    C->resetDOFManagers();

    // set that comms need to be refreshed
    // assumes that user will only call multiply again if they have changed
    // the values in A and or B
    // d_need_comms = true;
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::numericMultiplyReuse()
{
    PROFILE( "numericMultiplyReuse" );

    // start communication to build BRemote before doing anything
    if ( A->hasOffDiag() && d_need_comms ) {
        startBRemoteComm();
    }

    // First need to zero out coeffs in C so that remaining steps can be decoupled
    {
        lidx_t *C_rs = nullptr, *C_cols_loc = nullptr;
        gidx_t *C_cols                                 = nullptr;
        scalar_t *C_coeffs                             = nullptr;
        std::tie( C_rs, C_cols, C_cols_loc, C_coeffs ) = C_diag->getDataFields();
        auto tot_nnz                                   = C_diag->numberOfNonZeros();
        std::fill( C_coeffs, C_coeffs + tot_nnz, 0.0 );
        if ( !C_offd->isEmpty() ) {
            std::tie( C_rs, C_cols, C_cols_loc, C_coeffs ) = C_offd->getDataFields();
            tot_nnz                                        = C_offd->numberOfNonZeros();
            std::fill( C_coeffs, C_coeffs + tot_nnz, 0.0 );
        }
    }

    {
        PROFILE( "numericMultiplyReuse (A_diag)" );
        multiplyReuse( A_diag, B_diag, C_diag );
        multiplyReuse( A_diag, B_offd, C_offd );
    }

    if ( A->hasOffDiag() ) {
        PROFILE( "numericMultiplyReuse (A_offd)" );
        if ( d_need_comms ) {
            endBRemoteComm();
        }
        if ( BR_diag.get() != nullptr ) {
            multiplyReuse( A_offd, BR_diag, C_diag );
        }
        if ( BR_offd.get() != nullptr ) {
            multiplyReuse( A_offd, BR_offd, C_offd );
        }
    }

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
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    AMP_DEBUG_ASSERT( A_data != nullptr );
    AMP_DEBUG_ASSERT( B_data != nullptr );
    AMP_DEBUG_ASSERT( C_data != nullptr );

    if ( A_data->isEmpty() || B_data->isEmpty() ) {
        return;
    }

    const bool B_is_diag     = B_data->isDiag();
    constexpr bool use_dense = std::is_same_v<Accumulator, DenseAccumulator>;

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

    AMP_ASSERT( A_cols_loc != nullptr );
    if constexpr ( use_dense ) {
        AMP_ASSERT( B_cols_loc != nullptr ); // dense needs local cols
    }
    AMP_ASSERT( B_cols != nullptr || B_cols_loc != nullptr ); // otherwise just need one of them

    auto B_colmap          = B_data->getColumnMap();
    auto B_colmap_size     = B_data->numUniqueColumns();
    const auto B_first_col = B_data->beginCol();

    // DenseAcc's act on assembled blocks that have global columns removed
    // set up conversion for that case
    auto B_to_global =
        [B_cols, B_cols_loc, B_first_col, B_colmap, B_is_diag]( const lidx_t k ) -> gidx_t {
        if ( B_cols != nullptr ) {
            return B_cols[k];
        }
        return B_is_diag ? B_first_col + B_cols_loc[k] : B_colmap[B_cols_loc[k]];
    };

    // Create accumulator with appropriate capacity
    lidx_t acc_cap =
        use_dense ? ( B_is_diag ? B_data->numLocalColumns() : B_colmap_size ) : SPACC_SIZE;
    Accumulator acc( acc_cap, B_first_col );

    // Finally, after all the setup do the actual computation
    if constexpr ( IsSymbolic ) {
        PROFILE( use_dense ? "multiply (symbolic -- dense acc)" :
                             "multiply (symbolic -- sparse acc)" );
        // If this is a symbolic call just count NZ and write to
        // rs field in C
        for ( lidx_t row = 0; row < d_num_rows; ++row ) {
            // get rows in B block from the A_diag column indices
            for ( lidx_t j = A_rs[row]; j < A_rs[row + 1]; ++j ) {
                const auto Acl = A_cols_loc[j];
                // then row of C is union of those B row nz patterns
                for ( lidx_t k = B_rs[Acl]; k < B_rs[Acl + 1]; ++k ) {
                    const auto gbl = B_to_global( k );
                    acc.insert_or_append( gbl );
                }
            }
            // write out row length and clear accumulator
            C_rs[row] += acc.num_inserted;
            acc.clear();
        }
        C_data->setNNZ( true );
#if CSRSPGEMM_REPORT_SPACC_STATS
        if ( !use_dense && ( acc.total_collisions > 0 || acc.total_grows > 0 ) ) {
            AMP::pout << "\nSparseAcc stats:\n"
                      << "  Insertions: " << acc.total_inserted << "\n"
                      << "  Collisions: " << acc.total_collisions << "\n"
                      << "      Probes: " << acc.total_probe_steps << "\n"
                      << "      Clears: " << acc.total_clears << "\n"
                      << "       Grows: " << acc.total_grows << "\n"
                      << std::endl;
        }
#endif
    } else {
        PROFILE( use_dense ? "multiply (numeric -- dense acc)" :
                             "multiply (numeric -- sparse acc)" );
        // Otherwise, for numeric call write directly into C by
        // passing pointers into cols and coeffs fields as workspace
        // for the accumulator
        for ( lidx_t row = 0; row < d_num_rows; ++row ) {
            auto col_space = &C_cols[C_rs[row]];
            auto val_space = &C_coeffs[C_rs[row]];
            // get rows in B block from the A column indices
            for ( lidx_t j = A_rs[row]; j < A_rs[row + 1]; ++j ) {
                const auto Acl  = A_cols_loc[j];
                const auto Aval = A_coeffs[j];
                // then row of C is union of those B row nz patterns
                for ( lidx_t k = B_rs[Acl]; k < B_rs[Acl + 1]; ++k ) {
                    const auto gbl = B_to_global( k );
                    acc.insert_or_append( gbl, Aval * B_coeffs[k], col_space, val_space );
                }
            }
            // Clear accumulator to prepare for next row
            acc.clear();
        }
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
template<class Accumulator, bool IsDiag, bool IsSymbolic>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::multiplyFused(
    std::shared_ptr<DiagMatrixData> B_data,
    std::shared_ptr<DiagMatrixData> BR_data,
    std::shared_ptr<DiagMatrixData> C_data )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    if ( ( A_diag->isEmpty() || B_data->isEmpty() ) &&
         ( A_offd->isEmpty() || BR_data->isEmpty() ) ) {
        return;
    }

    constexpr bool use_dense = std::is_same_v<Accumulator, DenseAccumulator>;
    static_assert(use_dense == IsDiag);

    // all fields from blocks involved
    lidx_t *Ad_rs = nullptr, *Ad_cols_loc = nullptr;
    lidx_t *Ao_rs = nullptr, *Ao_cols_loc = nullptr;
    gidx_t *Ad_cols = nullptr, *Ao_cols = nullptr;
    scalar_t *Ad_coeffs = nullptr, *Ao_coeffs = nullptr;

    lidx_t *B_rs = nullptr, *B_cols_loc = nullptr;
    lidx_t *BR_rs = nullptr, *BR_cols_loc = nullptr;
    gidx_t *B_cols = nullptr, *BR_cols = nullptr;
    scalar_t *B_coeffs = nullptr, *BR_coeffs = nullptr;

    lidx_t *C_rs = nullptr, *C_cols_loc = nullptr;
    gidx_t *C_cols     = nullptr;
    scalar_t *C_coeffs = nullptr;

    // Extract available fields
    std::tie( Ad_rs, Ad_cols, Ad_cols_loc, Ad_coeffs ) = A_diag->getDataFields();
    std::tie( Ao_rs, Ao_cols, Ao_cols_loc, Ao_coeffs ) = A_offd->getDataFields();
    std::tie( B_rs, B_cols, B_cols_loc, B_coeffs )     = B_data->getDataFields();
    std::tie( BR_rs, BR_cols, BR_cols_loc, BR_coeffs ) = BR_data->getDataFields();
    std::tie( C_rs, C_cols, C_cols_loc, C_coeffs )     = C_data->getDataFields();

    AMP_ASSERT( Ad_cols_loc != nullptr && Ao_cols_loc != nullptr );
    AMP_ASSERT( B_data->isEmpty() || B_cols_loc != nullptr );
    AMP_ASSERT( BR_data->isEmpty() || BR_cols != nullptr );

    const auto first_col = C_data->beginCol();

    // The B blocks will have either local or global cols available
    // but generally not both. If only local available need conversion to global
    auto B_colmap    = B_offd->getColumnMap();
    auto B_to_global = [B_cols_loc, first_col, B_colmap]( const lidx_t k ) -> gidx_t {
        if constexpr ( IsDiag ) {
            return first_col + B_cols_loc[k];
        } else {
            return B_colmap[B_cols_loc[k]];
        }
    };

    // Create accumulator with appropriate capacity
    lidx_t acc_cap = use_dense ? B_data->numLocalColumns() : SPACC_SIZE;
    Accumulator acc( acc_cap, first_col );

    if constexpr ( IsSymbolic ) {
        PROFILE( use_dense ? "multiplyFusedSymbolic(dense acc)" :
                             "multiplyFusedSymbolic(sparse acc)" );
        // If this is a symbolic call just count NZ and write to
        // rs field in C
        for ( lidx_t row = 0; row < d_num_rows; ++row ) {
            // get rows in B block from the Ad column indices
            for ( lidx_t j = Ad_rs[row]; j < Ad_rs[row + 1]; ++j ) {
                const auto Acl = Ad_cols_loc[j];
                // then row of C is union of those B row nz patterns
                for ( lidx_t k = B_rs[Acl]; k < B_rs[Acl + 1]; ++k ) {
                    acc.insert_or_append( B_to_global( k ) );
                }
            }
            // get rows in BR block from the Ao column indices
            for ( lidx_t j = Ao_rs[row]; j < Ao_rs[row + 1]; ++j ) {
                const auto Acl = Ao_cols_loc[j];
                // then row of C is union of those B row nz patterns
                for ( lidx_t k = BR_rs[Acl]; k < BR_rs[Acl + 1]; ++k ) {
                    acc.insert_or_append( BR_cols[k] );
                }
            }
            // write out row length and clear accumulator
            C_rs[row] += acc.num_inserted;
            acc.clear();
        }
        C_data->setNNZ( true );
#if CSRSPGEMM_REPORT_SPACC_STATS
        if ( !use_dense && ( acc.total_collisions > 0 || acc.total_grows > 0 ) ) {
            AMP::pout << "\nSparseAcc stats:\n"
                      << "  Insertions: " << acc.total_inserted << "\n"
                      << "  Collisions: " << acc.total_collisions << "\n"
                      << "      Probes: " << acc.total_probe_steps << "\n"
                      << "      Clears: " << acc.total_clears << "\n"
                      << "       Grows: " << acc.total_grows << "\n"
                      << std::endl;
        }
#endif
    } else {
        PROFILE( use_dense ? "multiplyFusedNumeric(dense acc)" :
                             "multiplyFusedNumeric(sparse acc)" );
        // Otherwise, for numeric call write directly into C by
        // passing pointers into cols and coeffs fields as workspace
        // for the accumulator
        for ( lidx_t row = 0; row < d_num_rows; ++row ) {
            auto col_space = &C_cols[C_rs[row]];
            auto val_space = &C_coeffs[C_rs[row]];
            // get rows in B block from the Ad column indices
            for ( lidx_t j = Ad_rs[row]; j < Ad_rs[row + 1]; ++j ) {
                const auto Acl = Ad_cols_loc[j];
                // then row of C is union of those B row nz patterns
                for ( lidx_t k = B_rs[Acl]; k < B_rs[Acl + 1]; ++k ) {
                    const auto abVal = Ad_coeffs[j] * B_coeffs[k];
                    acc.insert_or_append( B_to_global( k ), abVal, col_space, val_space );
                }
            }
            // get rows in BR block from the Ao column indices
            for ( lidx_t j = Ao_rs[row]; j < Ao_rs[row + 1]; ++j ) {
                const auto Acl = Ao_cols_loc[j];
                // then row of C is union of those B row nz patterns
                for ( lidx_t k = BR_rs[Acl]; k < BR_rs[Acl + 1]; ++k ) {
                    const auto abVal = Ao_coeffs[j] * BR_coeffs[k];
                    acc.insert_or_append( BR_cols[k], abVal, col_space, val_space );
                }
            }
            // Clear accumulator to prepare for next row
            acc.clear();
        }
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::multiplyReuse(
    std::shared_ptr<DiagMatrixData> A_data,
    std::shared_ptr<DiagMatrixData> B_data,
    std::shared_ptr<DiagMatrixData> C_data )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    AMP_DEBUG_ASSERT( A_data != nullptr );
    AMP_DEBUG_ASSERT( B_data != nullptr );
    AMP_DEBUG_ASSERT( C_data != nullptr );

    if ( A_data->isEmpty() || B_data->isEmpty() ) {
        return;
    }

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

    const bool is_diag   = C_data->isDiag();
    auto B_colmap        = B_data->getColumnMap();
    auto C_colmap        = C_data->getColumnMap();
    const auto first_col = C_data->beginCol();

    AMP_ASSERT( A_cols_loc != nullptr );
    AMP_ASSERT( B_cols_loc != nullptr || B_cols != nullptr );
    AMP_ASSERT( C_cols_loc != nullptr );

    // columns already exist in C, so need way to convert them to position to write
    auto diag_search = [B_cols,
                        B_cols_loc]( lidx_t loc, lidx_t row_len, lidx_t *col_space ) -> lidx_t {
        for ( lidx_t pos = 0; pos < row_len; ++pos ) {
            if ( col_space[pos] == loc ) {
                return pos;
            }
        }
        AMP_ERROR( "multiplyReuse diag_search failed" );
        return -1;
    };
    auto offd_search = [C_colmap]( gidx_t gbl, lidx_t row_len, lidx_t *col_space ) -> lidx_t {
        for ( lidx_t pos = 0; pos < row_len; ++pos ) {
            if ( C_colmap[col_space[pos]] == gbl ) {
                return pos;
            }
        }
        AMP_ERROR( "multiplyReuse offd_search failed" );
        return -1;
    };

    // Finally, after all the setup do the actual computation
    {
        PROFILE( "multiplyReuse" );
        // Otherwise, for numeric call write directly into C by
        // passing pointers into cols and coeffs fields as workspace
        // for the accumulator
        for ( lidx_t row = 0; row < d_num_rows; ++row ) {
            const auto row_len = C_rs[row + 1] - C_rs[row];
            auto col_space     = &C_cols_loc[C_rs[row]];
            auto val_space     = &C_coeffs[C_rs[row]];
            // get rows in B block from the A column indices
            for ( lidx_t j = A_rs[row]; j < A_rs[row + 1]; ++j ) {
                const auto Acl = A_cols_loc[j];
                // then row of C is union of those B row nz patterns
                for ( lidx_t k = B_rs[Acl]; k < B_rs[Acl + 1]; ++k ) {
                    const auto idx =
                        is_diag ?
                            diag_search( B_cols != nullptr ?
                                             static_cast<lidx_t>( B_cols[k] - first_col ) :
                                             B_cols_loc[k],
                                         row_len,
                                         col_space ) :
                            offd_search( B_cols != nullptr ? B_cols[k] : B_colmap[B_cols_loc[k]],
                                         row_len,
                                         col_space );
                    val_space[idx] += A_coeffs[j] * B_coeffs[k];
                }
            }
        }
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
template<class Accumulator>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::mergeDiag()
{
    PROFILE( "mergeDiag" );

    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    const auto first_col = C_diag->beginCol();

    // handle special case where C_diag_offd is empty
    if ( C_diag_offd.get() == nullptr || C_diag_offd->isEmpty() ) {
        C_diag->swapDataFields( *C_diag_diag );
        return;
    }

    // pull out fields from blocks to merge and row pointers from C_diag
    lidx_t *C_dd_rs, *C_od_rs, *C_rs;
    lidx_t *C_dd_cols_loc, *C_od_cols_loc;
    gidx_t *C_dd_cols, *C_od_cols;
    scalar_t *C_dd_coeffs, *C_od_coeffs;

    std::tie( C_dd_rs, C_dd_cols, C_dd_cols_loc, C_dd_coeffs ) = C_diag_diag->getDataFields();
    std::tie( C_od_rs, C_od_cols, C_od_cols_loc, C_od_coeffs ) = C_offd_diag->getDataFields();
    C_rs                                                       = C_diag->getRowStarts();

    // Create allocator with space for C_diag operations
    const auto acc_cap = C_diag->numLocalColumns();
    Accumulator acc( acc_cap, first_col );

    // loop over all rows and count unique NZ positions in each
    for ( lidx_t row = 0; row < d_num_rows; ++row ) {
        // Add C_diag_diag row to accumulator
        for ( lidx_t j = C_dd_rs[row]; j < C_dd_rs[row + 1]; ++j ) {
            acc.insert_or_append( C_dd_cols[j] );
        }
        // Add C_diag_offd row to accumulator
        for ( lidx_t j = C_od_rs[row]; j < C_od_rs[row + 1]; ++j ) {
            acc.insert_or_append( C_od_cols[j] );
        }
        // write out row length and clear accumulator
        C_rs[row] += acc.num_inserted;
        acc.clear();
    }

    // allocate space in matrix
    C_diag->setNNZ( true );

    // pull result fields out
    lidx_t *C_cols_loc;
    gidx_t *C_cols;
    scalar_t *C_coeffs;
    std::tie( C_rs, C_cols, C_cols_loc, C_coeffs ) = C_diag->getDataFields();

    // loop over all rows again and write columns/coeffs into allocated matrix
    for ( lidx_t row = 0; row < d_num_rows; ++row ) {
        auto cols = &C_cols[C_rs[row]];
        auto vals = &C_coeffs[C_rs[row]];
        // Add C_diag_diag row to accumulator
        for ( lidx_t j = C_dd_rs[row]; j < C_dd_rs[row + 1]; ++j ) {
            acc.insert_or_append( C_dd_cols[j], C_dd_coeffs[j], cols, vals );
        }
        // Add C_diag_offd row to accumulator
        for ( lidx_t j = C_od_rs[row]; j < C_od_rs[row + 1]; ++j ) {
            acc.insert_or_append( C_od_cols[j], C_od_coeffs[j], cols, vals );
        }
        // clear accumulator
        acc.clear();
    }

    C_diag_diag.reset();
    C_offd_diag.reset();
}

template<typename Policy, class Allocator, class DiagMatrixData>
template<class Accumulator>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::mergeOffd()
{
    PROFILE( "mergeOffd" );

    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    // handle special case where either C_diag_offd or C_offd_offd is empty
    if ( C_diag_offd.get() == nullptr && C_offd_offd.get() == nullptr ) {
        return;
    }
    if ( C_offd_offd.get() == nullptr || C_offd_offd->isEmpty() ) {
        C_offd->swapDataFields( *C_diag_offd );
        return;
    }
    if ( C_diag_offd.get() == nullptr || C_diag_offd->isEmpty() ) {
        C_offd->swapDataFields( *C_offd_offd );
        return;
    }

    // pull out fields from blocks to merge and row pointers from C_offd
    lidx_t *C_do_rs, *C_oo_rs, *C_rs;
    lidx_t *C_do_cols_loc, *C_oo_cols_loc;
    gidx_t *C_do_cols, *C_oo_cols;
    scalar_t *C_do_coeffs, *C_oo_coeffs;

    std::tie( C_do_rs, C_do_cols, C_do_cols_loc, C_do_coeffs ) = C_diag_offd->getDataFields();
    std::tie( C_oo_rs, C_oo_cols, C_oo_cols_loc, C_oo_coeffs ) = C_offd_offd->getDataFields();
    C_rs                                                       = C_offd->getRowStarts();

    // Create allocator with space for C_offd operations
    Accumulator acc( SPACC_SIZE, 0 );

    // loop over all rows and count unique NZ positions in each
    for ( lidx_t row = 0; row < d_num_rows; ++row ) {
        // Add C_diag_offd row to accumulator
        for ( lidx_t j = C_do_rs[row]; j < C_do_rs[row + 1]; ++j ) {
            acc.insert_or_append( C_do_cols[j] );
        }
        // Add C_offd_offd row to accumulator
        for ( lidx_t j = C_oo_rs[row]; j < C_oo_rs[row + 1]; ++j ) {
            acc.insert_or_append( C_oo_cols[j] );
        }
        // write out row length and clear accumulator
        C_rs[row] += acc.num_inserted;
        acc.clear();
    }

    // allocate space in matrix
    C_offd->setNNZ( true );

    // report accumulator stats if useful
#if CSRSPGEMM_REPORT_SPACC_STATS
    if ( acc.total_collisions > 0 || acc.total_grows > 0 ) {
        AMP::pout << "\nSparseAcc stats:\n"
                  << "  Insertions: " << acc.total_inserted << "\n"
                  << "  Collisions: " << acc.total_collisions << "\n"
                  << "      Probes: " << acc.total_probe_steps << "\n"
                  << "      Clears: " << acc.total_clears << "\n"
                  << "       Grows: " << acc.total_grows << "\n"
                  << std::endl;
    }
#endif
    // pull result fields out
    lidx_t *C_cols_loc;
    gidx_t *C_cols;
    scalar_t *C_coeffs;
    std::tie( C_rs, C_cols, C_cols_loc, C_coeffs ) = C_offd->getDataFields();

    // loop over all rows again and write columns/coeffs into allocated matrix
    for ( lidx_t row = 0; row < d_num_rows; ++row ) {
        auto cols = &C_cols[C_rs[row]];
        auto vals = &C_coeffs[C_rs[row]];
        // Add C_diag_offd row to accumulator
        for ( lidx_t j = C_do_rs[row]; j < C_do_rs[row + 1]; ++j ) {
            acc.insert_or_append( C_do_cols[j], C_do_coeffs[j], cols, vals );
        }
        // Add C_offd_offd row to accumulator
        for ( lidx_t j = C_oo_rs[row]; j < C_oo_rs[row + 1]; ++j ) {
            acc.insert_or_append( C_oo_cols[j], C_oo_coeffs[j], cols, vals );
        }
        // clear accumulator
        acc.clear();
    }

    C_diag_offd.reset();
    C_offd_offd.reset();
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

    PROFILE( "setupBRemoteComm" );

    using lidx_t = typename Policy::lidx_t;

    auto comm_size = comm.getSize();

    // 1. Query comm list info and get offd colmap
    auto comm_list            = A->getRightCommList();
    auto rows_per_rank_recv   = comm_list->getReceiveSizes();
    auto rows_per_rank_send   = comm_list->getSendSizes();
    auto B_last_rows          = comm_list->getPartition();
    const auto A_col_map_size = A_offd->numUniqueColumns();
    const auto A_col_map      = A_offd->getColumnMap();

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
    PROFILE( "endBRemoteComm" );

    d_recv_matrices = d_csr_comm.recvMatrices( 0, 0, 0, B->numGlobalColumns() );
    // BRemotes do not need any particular parameters object internally
    BR_diag = CSRLocalMatrixData<Policy, Allocator>::ConcatVertical(
        nullptr, d_recv_matrices, B->beginCol(), B->endCol(), true );
    BR_offd = CSRLocalMatrixData<Policy, Allocator>::ConcatVertical(
        nullptr, d_recv_matrices, B->beginCol(), B->endCol(), false );

    // comms are done and BR_{diag,offd} filled, deallocate send/recv blocks
    d_send_matrices.clear();
    d_recv_matrices.clear();

    // set flag that recv'd matrices are valid
    d_need_comms = false;
}

template<typename Policy, class Allocator, class DiagMatrixData>
bool CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::DenseAccumulator::contains(
    typename Policy::gidx_t gbl ) const
{
    using lidx_t = typename Policy::lidx_t;

    const auto loc = static_cast<lidx_t>( gbl - offset );
    const auto k   = flags[loc];
    return k != -1;
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::DenseAccumulator::
    insert_or_append( typename Policy::gidx_t gbl )
{
    using lidx_t = typename Policy::lidx_t;

    const auto loc = static_cast<lidx_t>( gbl - offset );
    const auto k   = flags[loc];
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
    insert_or_append( typename Policy::gidx_t gbl,
                      typename Policy::scalar_t val,
                      typename Policy::gidx_t *col_space,
                      typename Policy::scalar_t *val_space )
{
    using lidx_t = typename Policy::lidx_t;

    const auto loc = static_cast<lidx_t>( gbl - offset );
    const auto k   = flags[loc];
    if ( k == -1 ) {
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
uint16_t CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::SparseAccumulator::hash(
    typename Policy::gidx_t gbl ) const
{
    const uint16_t c0 = ( 506999 * gbl ) & 0xFFFF;
    const uint16_t c1 = ( gbl >> 16 ) & 0xFFFF;
    return ( c0 ^ c1 ) % capacity;
}

template<typename Policy, class Allocator, class DiagMatrixData>
bool CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::SparseAccumulator::contains(
    typename Policy::gidx_t gbl ) const
{
    auto pos = hash( gbl ), flag = flags[pos];
    if ( flag == 0xFFFF ) {
        // Location is empty, certainly not present
        return false;
    } else {
        // location occupied, linear probe to empty or gbl found
        do {
            if ( cols[flag] == gbl ) {
                return true;
            }
            pos  = ( pos + 1 ) % capacity;
            flag = flags[pos];
        } while ( flag != 0xFFFF );
    }
    // gbl never found, is not contained
    return false;
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::SparseAccumulator::
    insert_or_append( typename Policy::gidx_t gbl )
{
    if ( num_inserted == capacity ) {
        grow( cols.data() );
    }

    auto pos = hash( gbl ), flag = flags[pos];
    if ( flag == 0xFFFF ) {
        // Location is empty, append
        flags[pos] = num_inserted;
        if ( cols.size() <= num_inserted ) {
            cols.push_back( gbl );
        } else {
            cols[num_inserted] = gbl;
        }
        ++num_inserted;
#if CSRSPGEMM_REPORT_SPACC_STATS
        ++total_inserted;
#endif
    } else {
        // location occupied, linear probe to empty or gbl found
#if CSRSPGEMM_REPORT_SPACC_STATS
        if ( cols[flag] != gbl ) {
            ++total_collisions;
        }
#endif
        do {
            if ( cols[flag] == gbl ) {
                return;
            }
            pos  = ( pos + 1 ) % capacity;
            flag = flags[pos];
#if CSRSPGEMM_REPORT_SPACC_STATS
            ++total_probe_steps;
#endif
        } while ( flag != 0xFFFF );
        // gbl never found, have empty slot
        flags[pos] = num_inserted;
        if ( cols.size() <= num_inserted ) {
            cols.push_back( gbl );
        } else {
            cols[num_inserted] = gbl;
        }
        ++num_inserted;
#if CSRSPGEMM_REPORT_SPACC_STATS
        ++total_inserted;
#endif
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::SparseAccumulator::
    insert_or_append( typename Policy::gidx_t gbl,
                      typename Policy::scalar_t val,
                      typename Policy::gidx_t *col_space,
                      typename Policy::scalar_t *val_space )
{
    if ( num_inserted == capacity ) {
        grow( col_space );
    }

    auto pos = hash( gbl ), flag = flags[pos];
    if ( flag == 0xFFFF ) {
        // Location is empty, append
        flags[pos]              = num_inserted;
        col_space[num_inserted] = gbl;
        val_space[num_inserted] = val;
        ++num_inserted;
    } else {
        // location occupied, linear probe to empty or gbl found
        do {
            if ( col_space[flag] == gbl ) {
                val_space[flag] += val;
                return;
            }
            pos  = ( pos + 1 ) % capacity;
            flag = flags[pos];
        } while ( flag != 0xFFFF );
        // gbl never found, have empty slot
        flags[pos]              = num_inserted;
        col_space[num_inserted] = gbl;
        val_space[num_inserted] = val;
        ++num_inserted;
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::SparseAccumulator::grow(
    typename Policy::gidx_t *col_space )
{
#if CSRSPGEMM_REPORT_SPACC_STATS
    ++total_grows;
#endif

    uint16_t old_capacity = capacity;
    capacity *= 2;
    std::vector<uint16_t> old_flags( capacity, 0xFFFF );
    flags.swap( old_flags );
    for ( uint16_t n = 0; n < old_capacity; ++n ) {
        // insert all existing global local pairs into new flags vector
        // do it inline since cols/col_space should not be touched,
        // and nor should the num_inserted value
        const auto loc = old_flags[n];
        const auto gbl = col_space[loc];
        auto pos = hash( gbl ), flag = flags[pos];
        if ( flag == 0xFFFF ) {
            // Location is empty, set and carry on
            flags[pos] = loc;
        } else {
            // location occupied, linear probe to empty
            // Note: must find empty since grow is guaranteed to work with set
            // of unique columns
            do {
                pos  = ( pos + 1 ) % capacity;
                flag = flags[pos];
            } while ( flag != 0xFFFF );
            flags[pos] = loc;
        }
    }
}

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixSpGEMMHelperDefault<Policy, Allocator, DiagMatrixData>::SparseAccumulator::clear()
{
#if CSRSPGEMM_REPORT_SPACC_STATS
    total_clears++;
#endif
    num_inserted = 0;
    std::fill( flags.begin(), flags.end(), 0xFFFF );
}

} // namespace AMP::LinearAlgebra

#ifdef CSRSPGEMM_REPORT_SPACC_STATS
    #undef CSRSPGEMM_REPORT_SPACC_STATS
#endif
