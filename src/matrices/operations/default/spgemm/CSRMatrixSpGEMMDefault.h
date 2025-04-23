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
    template<bool SYMBOLIC>
    void multiplyLocal( std::shared_ptr<DiagMatrixData> B_data,
                        std::shared_ptr<DiagMatrixData> C_data,
                        lidx_t *nnz );
    void symbolicMultiplyLocal( std::shared_ptr<DiagMatrixData> B_data, std::vector<lidx_t> &nnz );
    void numericMultiplyLocal( std::shared_ptr<DiagMatrixData> B_data,
                               std::shared_ptr<DiagMatrixData> C_data );

    void setupBRemoteComm();
    void startBRemoteComm();
    void endBRemoteComm();

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
    // This is a single block for all columns in remote rows. It is much easier
    // to split these between C_diag and C_offd after receiving them
    std::shared_ptr<DiagMatrixData> BRemote;

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

    // Internal row accumlator class
    struct DenseAccumulator {
        DenseAccumulator( int capacity_ )
            : capacity( capacity_ ), num_inserted( 0 ), flags( capacity, -1 )
        {
        }

        void insert_or_append( lidx_t loc, gidx_t gbl )
        {
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

        void insert_or_append(
            lidx_t loc, gidx_t gbl, scalar_t val, gidx_t *col_space, scalar_t *val_space )
        {
            const auto k = flags[loc];
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

        void clear()
        {
            for ( int n = 0; n < num_inserted; ++n ) {
                flags[flag_inv[n]] = -1;
            }
            num_inserted = 0;
        }

        const lidx_t capacity;
        lidx_t num_inserted;
        std::vector<lidx_t> flags;
        std::vector<lidx_t> flag_inv;
        std::vector<gidx_t> cols;
    };
};

} // namespace AMP::LinearAlgebra

#endif
