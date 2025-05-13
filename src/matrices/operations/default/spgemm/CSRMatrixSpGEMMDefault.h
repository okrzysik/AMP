#ifndef included_AMP_CSRMatrixSpGEMMDefault
#define included_AMP_CSRMatrixSpGEMMDefault

#include "AMP/matrices/data/CSRMatrixCommunicator.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/utils/AMP_MPI.h"

#include <map>
#include <memory>
#include <vector>

namespace AMP::LinearAlgebra {

template<typename Policy, class Allocator, class DiagMatrixData>
class CSRMatrixSpGEMMHelperDefault
{
    using CSRData  = CSRMatrixData<Policy, Allocator, DiagMatrixData>;
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

public:
    CSRMatrixSpGEMMHelperDefault() = default;
    CSRMatrixSpGEMMHelperDefault( CSRData *A_, CSRData *B_, CSRData *C_, bool overlap_comms_ )
        : A( A_ ),
          B( B_ ),
          C( C_ ),
          A_diag( A->getDiagMatrix() ),
          A_offd( A->getOffdMatrix() ),
          B_diag( B->getDiagMatrix() ),
          B_offd( B->getOffdMatrix() ),
          C_diag( C->getDiagMatrix() ),
          C_offd( C->getOffdMatrix() ),
          d_overlap_comms( ( A->getComm().getSize() > 1 ) && overlap_comms_ ),
          comm( A->getComm() ),
          d_csr_comm( A->getRightCommList() ),
          d_need_comms( true )
    {
        AMP_DEBUG_INSIST(
            comm == B->getComm() && comm == C->getComm(),
            "CSRMatrixSpGEMMHelperDefault: All three matrices must have the same communicator" );
    }

    ~CSRMatrixSpGEMMHelperDefault() = default;

    void symbolicMultiply();
    void numericMultiply();
    void numericMultiplyReuse();

protected:
    void symbolicMultiply_NonOverlapped();
    void numericMultiply_NonOverlapped();
    void symbolicMultiply_Overlapped();
    void numericMultiply_Overlapped();

    template<class Accumulator, bool IsSymbolic>
    void multiply( std::shared_ptr<DiagMatrixData> A_data,
                   std::shared_ptr<DiagMatrixData> B_data,
                   std::shared_ptr<DiagMatrixData> C_data,
                   lidx_t *ctr );

    template<class Accumulator, bool IsSymbolic>
    void multiplyFused( std::shared_ptr<DiagMatrixData> B0_data,
                        std::shared_ptr<DiagMatrixData> B1_data,
                        std::shared_ptr<DiagMatrixData> C_data );

    void setupBRemoteComm();
    void startBRemoteComm();
    void endBRemoteComm();

    template<class Accumulator>
    void mergeDiag();

    template<class Accumulator>
    void mergeOffd();

    // Useful constants for supporting operations
    // Fill factor used to estimate size of output matrix NZs prior to
    // symbolic phase. This is applied *before* BRemote is gathered
    // from other ranks, so make it a little higher than perhaps expected
    static constexpr scalar_t C_FILL_MULT   = 1.5;
    static constexpr scalar_t C_FILL_ADD    = 0.75;
    static constexpr scalar_t C_GROW_FACTOR = 1.5;
    // default starting size for sparse accumulators
    static constexpr lidx_t SPACC_SIZE = 256;

    // Matrix data of operands and output
    // these are non-owning pointers
    CSRData *A;
    CSRData *B;
    CSRData *C;

    // diag and offd blocks of input matrices
    std::shared_ptr<DiagMatrixData> A_diag;
    std::shared_ptr<DiagMatrixData> A_offd;
    std::shared_ptr<DiagMatrixData> B_diag;
    std::shared_ptr<DiagMatrixData> B_offd;

    // Matrix data formed from remote rows of B that get pulled to each process
    std::shared_ptr<DiagMatrixData> BR_diag;
    std::shared_ptr<DiagMatrixData> BR_offd;

    // Blocks of C matrix
    std::shared_ptr<DiagMatrixData> C_diag;
    std::shared_ptr<DiagMatrixData> C_offd;

    // flag for whether overlapped communication/computation should be done
    bool d_overlap_comms;

    // Communicator
    AMP_MPI comm;
    CSRMatrixCommunicator<Policy, Allocator, DiagMatrixData> d_csr_comm;
    bool d_need_comms;

    // To overlap comms and calcs it is easiest to form the output in four
    // blocks and merge them together at the end
    std::shared_ptr<DiagMatrixData> C_diag_diag; // from A_diag * B_diag
    std::shared_ptr<DiagMatrixData> C_diag_offd; // from A_diag * B_offd
    std::shared_ptr<DiagMatrixData> C_offd_diag; // from A_offd * BR_diag
    std::shared_ptr<DiagMatrixData> C_offd_offd; // from A_offd * BR_offd

    // The following all support the communication needed to build BRemote
    // these are worth preserving to allow repeated SpGEMMs to re-use the
    // structure with potentially new coefficients

    // struct to hold fields that are needed in both
    // the "source" and "destination" perspectives
    struct SpGEMMCommInfo {
        SpGEMMCommInfo() = default;
        SpGEMMCommInfo( int numrow_ ) : numrow( numrow_ ) {}
        // number of rows to send or receive
        int numrow;
        // ids of rows to send/receive
        std::vector<gidx_t> rowids;
        // number of non-zeros in those rows
        std::vector<lidx_t> rownnz;
    };

    // Source information, things expected from other ranks
    std::map<int, SpGEMMCommInfo> d_src_info;

    // Destination information, things sent to other ranks
    std::map<int, SpGEMMCommInfo> d_dest_info;

    std::map<int, std::shared_ptr<DiagMatrixData>> d_send_matrices;
    std::map<int, std::shared_ptr<DiagMatrixData>> d_recv_matrices;

    // Internal row accumlator classes
    struct DenseAccumulator {
        DenseAccumulator( int capacity_, gidx_t offset_ )
            : capacity( capacity_ ),
              offset( offset_ ),
              num_inserted( 0 ),
              total_inserted( 0 ),
              total_collisions( 0 ),
              total_probe_steps( 0 ),
              total_clears( 0 ),
              total_grows( 0 ),
              flags( capacity, -1 )
        {
        }

        void insert_or_append( gidx_t gbl );
        void insert_or_append(
            gidx_t gbl, scalar_t val, gidx_t *col_space, scalar_t *val_space, lidx_t max_pos );
        void clear();

        const lidx_t capacity;
        const gidx_t offset;
        lidx_t num_inserted;
        size_t total_inserted;
        size_t total_collisions;
        size_t total_probe_steps;
        size_t total_clears;
        size_t total_grows;
        std::vector<lidx_t> flags;
        std::vector<lidx_t> flag_inv;
        std::vector<gidx_t> cols;
    };

    struct SparseAccumulator {
        SparseAccumulator( int capacity_, gidx_t offset_ )
            : capacity( capacity_ ),
              offset( offset_ ),
              num_inserted( 0 ),
              total_inserted( 0 ),
              total_collisions( 0 ),
              total_probe_steps( 0 ),
              total_clears( 0 ),
              total_grows( 0 ),
              flags( capacity, 0xFFFF )
        {
            AMP_DEBUG_ASSERT( capacity > 1 );
        }

        uint16_t hash( gidx_t gbl ) const;
        void insert_or_append( gidx_t gbl );
        void insert_or_append(
            gidx_t gbl, scalar_t val, gidx_t *col_space, scalar_t *val_space, lidx_t max_pos );
        void clear();

        uint16_t capacity;
        const gidx_t offset;
        uint16_t num_inserted;
        size_t total_inserted;
        size_t total_collisions;
        size_t total_probe_steps;
        size_t total_clears;
        size_t total_grows;
        std::vector<uint16_t> flags;
        std::vector<gidx_t> cols;

    private:
        void grow( gidx_t *col_space );
    };
};

} // namespace AMP::LinearAlgebra

#endif
