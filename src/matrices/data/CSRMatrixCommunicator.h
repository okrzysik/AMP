#ifndef included_AMP_CSRMatrixCommunicator_h
#define included_AMP_CSRMatrixCommunicator_h

#include "AMP/matrices/data/CSRLocalMatrixData.h"
#include "AMP/utils/AMP_MPI.h"

#include <map>
#include <vector>

namespace AMP::LinearAlgebra {

template<typename Policy,
         class Allocator      = AMP::HostAllocator<void>,
         class DiagMatrixData = CSRLocalMatrixData<Policy, Allocator>>
class CSRMatrixCommunicator
{
public:
    using gidx_t   = typename Policy::gidx_t;
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;
    static_assert( std::is_same_v<typename Allocator::value_type, void> );
    using gidxAllocator_t =
        typename std::allocator_traits<Allocator>::template rebind_alloc<gidx_t>;
    using lidxAllocator_t =
        typename std::allocator_traits<Allocator>::template rebind_alloc<lidx_t>;
    using scalarAllocator_t =
        typename std::allocator_traits<Allocator>::template rebind_alloc<scalar_t>;

    CSRMatrixCommunicator() = default;
    CSRMatrixCommunicator( AMP_MPI &comm )
        : d_comm( comm ), d_send_called( false ), d_num_sources( 0 )
    {
    }

    void sendMatrices( std::map<int, std::shared_ptr<DiagMatrixData>> matrices );
    std::map<int, std::shared_ptr<DiagMatrixData>>
    recvMatrices( gidx_t first_row, gidx_t last_row, gidx_t first_col, gidx_t last_col );

protected:
    AMP_MPI d_comm;
    bool d_send_called;
    int d_num_sources;

    std::vector<AMP_MPI::Request> d_send_requests;
    std::vector<AMP_MPI::Request> d_recv_requests;

    // tags for each type of message to send/recv
    static constexpr int ROW_TAG   = 5600;
    static constexpr int COL_TAG   = 5601;
    static constexpr int COEFF_TAG = 5602;
};
}

#endif
