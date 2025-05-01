#ifndef included_AMP_CSRMatrixSpGEMMDefault
#define included_AMP_CSRMatrixSpGEMMDefault

#include "AMP/matrices/data/CSRMatrixCommunicator.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/utils/AMP_MPI.h"

#include <map>
#include <memory>
#include <set>
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
    CSRMatrixSpGEMMHelperDefault( CSRData *A_, CSRData *B_, CSRData *C_ )
        : A( A_ ),
          B( B_ ),
          C( C_ ),
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

protected:
    template<class Accumulator, bool IsSymbolic>
    void multiply( std::shared_ptr<DiagMatrixData> A_data,
                   std::shared_ptr<DiagMatrixData> B_data,
                   std::shared_ptr<DiagMatrixData> C_data );

    void setupBRemoteComm();
    void startBRemoteComm();
    void endBRemoteComm();

    template<class Accumulator>
    void mergeDiag();

    template<class Accumulator>
    void mergeOffd();

    // This only (re-)fills the coefficients in BRemote
    // the comm info and symbolic creation must have already happened
    // void fillBRemoteNumeric();

    // Matrix data of operands and output
    // these are non-owning pointers
    CSRData *A;
    CSRData *B;
    CSRData *C;

    // Communicator
    AMP_MPI comm;
    CSRMatrixCommunicator<Policy, Allocator, DiagMatrixData> d_csr_comm;
    bool d_need_comms;

    // Matrix data formed from remote rows of B that get pulled to each process
    std::shared_ptr<DiagMatrixData> BR_diag;
    std::shared_ptr<DiagMatrixData> BR_offd;

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
        DenseAccumulator( int capacity_ )
            : capacity( capacity_ ), num_inserted( 0 ), flags( capacity, -1 )
        {
        }

        void insert_or_append( lidx_t loc, gidx_t gbl );
        void insert_or_append( lidx_t loc,
                               gidx_t gbl,
                               scalar_t val,
                               gidx_t *col_space,
                               scalar_t *val_space,
                               lidx_t max_pos );
        void clear();

        const lidx_t capacity;
        lidx_t num_inserted;
        std::vector<lidx_t> flags;
        std::vector<lidx_t> flag_inv;
        std::vector<gidx_t> cols;
    };

    struct SparseAccumulator {
        SparseAccumulator( int capacity_ ) : capacity( capacity_ ), num_inserted( 0 ) {}

        void insert_or_append( lidx_t, gidx_t gbl );
        void insert_or_append( lidx_t,
                               gidx_t gbl,
                               scalar_t val,
                               gidx_t *col_space,
                               scalar_t *val_space,
                               lidx_t max_pos );
        void clear();

        const lidx_t capacity;
        lidx_t num_inserted;
        std::unordered_map<gidx_t, lidx_t> kv;
    };
};

} // namespace AMP::LinearAlgebra

#endif
