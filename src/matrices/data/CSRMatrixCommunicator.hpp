#ifndef included_AMP_CSRMatrixCommunicator_hpp
#define included_AMP_CSRMatrixCommunicator_hpp

#include "AMP/matrices/data/CSRMatrixCommunicator.h"

namespace AMP::LinearAlgebra {

template<typename Policy, class Allocator, class DiagMatrixData>
void CSRMatrixCommunicator<Policy, Allocator, DiagMatrixData>::sendMatrices(
    std::map<int, std::shared_ptr<DiagMatrixData>> matrices )
{
    AMP_ASSERT( d_send_requests.size() == 0 );
    std::vector<int> all_dests( d_comm.getSize(), 0 );
    for ( auto it : matrices ) {
        const int dest     = it.first;
        all_dests[dest]    = 1;
        auto matrix        = it.second;
        const auto num_rs  = matrix->d_num_rows + 1;
        const auto num_nnz = matrix->d_nnz;
        d_send_requests.emplace_back(
            d_comm.Isend( matrix->d_row_starts.get(), num_rs, dest, ROW_TAG ) );
        d_send_requests.emplace_back(
            d_comm.Isend( matrix->d_cols.get(), num_nnz, dest, COL_TAG ) );
        d_send_requests.emplace_back(
            d_comm.Isend( matrix->d_coeffs.get(), num_nnz, dest, COEFF_TAG ) );
    }
    auto all_sources = d_comm.allToAll( all_dests );
    d_num_sources    = std::reduce( all_sources.begin(), all_sources.end() );
    d_send_called    = true;
}

template<typename Policy, class Allocator, class DiagMatrixData>
std::map<int, std::shared_ptr<DiagMatrixData>>
CSRMatrixCommunicator<Policy, Allocator, DiagMatrixData>::recvMatrices(
    typename Policy::gidx_t first_row,
    typename Policy::gidx_t last_row,
    typename Policy::gidx_t first_col,
    typename Policy::gidx_t last_col )
{
    using lidx_t = typename Policy::lidx_t;
    using gidx_t = typename Policy::gidx_t;

    AMP_INSIST( d_send_called,
                "CSRMatrixCommunicator::sendMatrices must be called before recvMatrices" );

    std::map<int, std::shared_ptr<DiagMatrixData>> blocks;
    const auto mem_loc = AMP::Utilities::getAllocatorMemoryType<Allocator>();

    // there are d_num_sources matrices to recieve
    // always sent in order row_starts, cols, coeffs
    // start with probe on any source with ROW_TAG
    for ( int ns = 0; ns < d_num_sources; ++ns ) {
        auto [source, tag, num_bytes] = d_comm.probe( -1, ROW_TAG );
        AMP_ASSERT( tag == ROW_TAG );
        // remember row_starts has extra entry
        const lidx_t num_rows = ( num_bytes / sizeof( lidx_t ) ) - 1;
        // if last_row is zero then choose based on num_rows,
        // otherwise test that recv'd matrix is valid with layout
        gidx_t fr, lr;
        if ( last_row == 0 ) {
            fr = 0;
            lr = static_cast<gidx_t>( num_rows );
        } else {
            AMP_INSIST( num_rows == static_cast<lidx_t>( last_row - first_row ),
                        "Received matrix with invalid layout" );
            fr = first_row;
            lr = last_row;
        }
        auto [it, inserted] =
            blocks.insert( { source,
                             std::make_shared<DiagMatrixData>(
                                 nullptr, mem_loc, fr, lr, first_col, last_col, false ) } );
        AMP_ASSERT( inserted );
        auto block = ( *it ).second;
        // matrix now exists and has row_starts buffer, recv it and trigger allocations
        d_comm.recv( block->d_row_starts.get(), num_rows + 1, source, ROW_TAG );
        block->setNNZ( false );
        // buffers for cols and coeffs now allocated, recv them and continue to next probe
        d_comm.recv( block->d_cols.get(), block->d_nnz, source, COL_TAG );
        d_comm.recv( block->d_coeffs.get(), block->d_nnz, source, COEFF_TAG );
    }

    // enaure that any outstanding sends complete
    d_comm.waitAll( static_cast<int>( d_send_requests.size() ), d_send_requests.data() );

    // comm done, reset send flag in case this gets re-used
    d_send_called = false;

    return blocks;
}

} // namespace AMP::LinearAlgebra

#endif
