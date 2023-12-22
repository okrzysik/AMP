#include "AMP/matrices/data/hypre/HYPREMatrixAdaptor.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/utils/AMP_MPI.h"

namespace AMP::LinearAlgebra {

HypreMatrixAdaptor::HypreMatrixAdaptor( std::shared_ptr<MatrixData> matrixData )
{
    HYPRE_BigInt first_row = static_cast<HYPRE_BigInt>( matrixData->beginRow() );
    HYPRE_BigInt last_row  = static_cast<HYPRE_BigInt>( matrixData->endRow() - 1 );
    auto comm              = matrixData->getComm().getCommunicator();

    HYPRE_IJMatrixCreate( comm, first_row, last_row, first_row, last_row, &d_matrix );
    HYPRE_IJMatrixSetObjectType( d_matrix, HYPRE_PARCSR );

    HYPRE_Int *nnz_per_row = nullptr;
    HYPRE_BigInt *csr_ja   = nullptr;
    HYPRE_Real *csr_aa     = nullptr;

    auto csrData = std::dynamic_pointer_cast<CSRMatrixData>( matrixData );
    if ( csrData ) {
        AMP_ERROR( "Not implemented" );
    } else {
        AMP_ERROR( "Not implemented" );
    }

    AMP_INSIST( nnz_per_row && csr_ja && csr_aa, "nnz_per_row, csr_ja, csr_aa cannot be NULL" );
    initializeHypreMatrix( first_row, last_row, nnz_per_row, csr_ja, csr_aa );
}

HypreMatrixAdaptor::~HypreMatrixAdaptor() { HYPRE_IJMatrixDestroy( d_matrix ); }

static void
set_row_ids_( HYPRE_BigInt const first_row, HYPRE_BigInt const nrows, HYPRE_BigInt *row_ids )
{
    AMP_ERROR( "Not implemented" );
}

void HypreMatrixAdaptor::initializeHypreMatrix( HYPRE_BigInt first_row,
                                                HYPRE_BigInt last_row,
                                                HYPRE_Int *const nnz_per_row,
                                                HYPRE_BigInt *const csr_ja,
                                                HYPRE_Real *const csr_aa )
{
    const auto nrows = last_row - first_row + 1;

    HYPRE_IJMatrixSetRowSizes( d_matrix, nnz_per_row );
    HYPRE_IJMatrixSetMaxOffProcElmts( d_matrix, 0 );

    // The next 2 lines affect efficiency and should be resurrected at some point
    //  set_row_location_(d_first_row, d_last_row, nrows, nnz_per_row, csr_ia, csr_ja,
    //  number_of_local_cols, number_of_remote_cols ); HYPRE_IJMatrixSetDiagOffdSizes( d_matrix,
    //  number_of_local_cols.data(), number_of_remote_cols.data() );

    HYPRE_IJMatrixInitialize( d_matrix );

    HYPRE_BigInt *row_ids = nullptr;
    AMP_INSIST( row_ids, "row_ids cannot be NULL" );
    set_row_ids_( first_row, nrows, row_ids );
    HYPRE_IJMatrixSetValues( d_matrix, nrows, nnz_per_row, row_ids, csr_ja, csr_aa );
    HYPRE_IJMatrixAssemble( d_matrix );
}

} // namespace AMP::LinearAlgebra
