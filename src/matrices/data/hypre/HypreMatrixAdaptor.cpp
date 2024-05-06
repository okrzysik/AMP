#include "AMP/matrices/data/hypre/HypreMatrixAdaptor.h"
#include "AMP/AMP_TPLs.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/data/hypre/HypreCSRPolicy.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"

#ifdef USE_CUDA
    #include "AMP/utils/cuda/CudaAllocator.h"
#endif

#include <numeric>

namespace AMP::LinearAlgebra {

HypreMatrixAdaptor::HypreMatrixAdaptor( std::shared_ptr<MatrixData> matrixData )
{
    int ierr;
    char hypre_mesg[100];

    HYPRE_BigInt firstRow = static_cast<HYPRE_BigInt>( matrixData->beginRow() );
    HYPRE_BigInt lastRow  = static_cast<HYPRE_BigInt>( matrixData->endRow() - 1 );
    auto comm             = matrixData->getComm().getCommunicator();

    HYPRE_IJMatrixCreate( comm, firstRow, lastRow, firstRow, lastRow, &d_matrix );
    HYPRE_IJMatrixSetObjectType( d_matrix, HYPRE_PARCSR );

    HYPRE_IJMatrixSetMaxOffProcElmts( d_matrix, 0 );

    auto csrData = std::dynamic_pointer_cast<CSRMatrixData<HypreCSRPolicy>>( matrixData );
    if ( csrData ) {

        auto [nnz_d, cols_d, cols_loc_d, coeffs_d]     = csrData->getCSRDiagData();
        auto [nnz_od, cols_od, cols_loc_od, coeffs_od] = csrData->getCSROffDiagData();

        auto nnz_per_row_d  = static_cast<HYPRE_Int *>( nnz_d );
        const auto csr_ja_d = static_cast<HYPRE_BigInt const *>( cols_d );
        const auto csr_aa_d = static_cast<HYPRE_Real const *>( coeffs_d );

        auto nnz_per_row_od  = static_cast<HYPRE_Int *>( nnz_od );
        const auto csr_ja_od = static_cast<HYPRE_BigInt const *>( cols_od );
        const auto csr_aa_od = static_cast<HYPRE_Real const *>( coeffs_od );

        AMP_INSIST( nnz_d && csr_ja_d && csr_aa_d, "diagonal block layout cannot be NULL" );
        initializeHypreMatrix( firstRow,
                               lastRow,
                               csrData->hasOffDiag(),
                               nnz_per_row_d,
                               csr_ja_d,
                               csr_aa_d,
                               nnz_per_row_od,
                               csr_ja_od,
                               csr_aa_od );

    } else {

        // figure out how to incorporate this
        // HYPRE_IJMatrixSetRowSizes( d_matrix, nnz_per_row );

        HYPRE_IJMatrixInitialize( d_matrix );

        // iterate over all rows
        for ( auto i = firstRow; i <= lastRow; ++i ) {

            std::vector<size_t> cols;
            std::vector<double> values;

            matrixData->getRowByGlobalID( i, cols, values );
            std::vector<HYPRE_BigInt> hypre_cols( cols.size() );
            std::copy( cols.begin(), cols.end(), hypre_cols.begin() );

            const int nrows  = 1;
            const auto irow  = i;
            const auto ncols = cols.size();

            ierr = HYPRE_IJMatrixSetValues( d_matrix,
                                            nrows,
                                            (HYPRE_Int *) &ncols,
                                            (HYPRE_BigInt *) &irow,
                                            hypre_cols.data(),
                                            (const HYPRE_Real *) values.data() );
            HYPRE_DescribeError( ierr, hypre_mesg );
        }

        HYPRE_IJMatrixAssemble( d_matrix );
    }
}

HypreMatrixAdaptor::~HypreMatrixAdaptor() { HYPRE_IJMatrixDestroy( d_matrix ); }

static void
set_row_ids_( HYPRE_BigInt const first_row, HYPRE_BigInt const nrows, HYPRE_BigInt *row_ids )
{
    AMP_ASSERT( row_ids );
    std::iota( row_ids, row_ids + nrows, first_row );
}

void HypreMatrixAdaptor::initializeHypreMatrix( HYPRE_BigInt first_row,
                                                HYPRE_BigInt last_row,
                                                bool has_off_diag,
                                                HYPRE_Int *const nnz_per_row_d,
                                                HYPRE_BigInt const *const csr_ja_d,
                                                HYPRE_Real const *const csr_aa_d,
                                                HYPRE_Int *const nnz_per_row_od,
                                                HYPRE_BigInt const *const csr_ja_od,
                                                HYPRE_Real const *const csr_aa_od )
{
    const auto nrows = last_row - first_row + 1;

#ifdef USE_CUDA
    AMP::CudaManagedAllocator<HYPRE_BigInt> managedAllocator;
#endif

    HYPRE_IJMatrixSetDiagOffdSizes( d_matrix, nnz_per_row_d, nnz_per_row_od );

    HYPRE_IJMatrixInitialize( d_matrix );

    auto memType = AMP::Utilities::getMemoryType( csr_ja_d );

    HYPRE_BigInt *row_ids_p = nullptr;
    std::vector<HYPRE_BigInt> row_ids;
    if ( memType <= AMP::Utilities::MemoryType::host ) {
        row_ids.resize( nrows );
        row_ids_p = row_ids.data();
    } else if ( memType == AMP::Utilities::MemoryType::managed ) {

#ifdef USE_CUDA
        row_ids_p = managedAllocator.allocate( nrows );
#endif

    } else if ( memType == AMP::Utilities::MemoryType::device ) {
        AMP_ERROR( "Currently only implemented for host accessible memory" );
    }

    set_row_ids_( first_row, nrows, row_ids_p );
    HYPRE_IJMatrixSetValues( d_matrix, nrows, nnz_per_row_d, row_ids_p, csr_ja_d, csr_aa_d );
    if ( has_off_diag ) {
        HYPRE_IJMatrixSetValues( d_matrix, nrows, nnz_per_row_od, row_ids_p, csr_ja_od, csr_aa_od );
    }
    HYPRE_IJMatrixAssemble( d_matrix );

    if ( memType == AMP::Utilities::MemoryType::managed ) {
#ifdef USE_CUDA
        managedAllocator.deallocate( row_ids_p, nrows );
#endif
    }
}

} // namespace AMP::LinearAlgebra
