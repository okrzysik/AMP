#include "AMP/matrices/data/hypre/HypreMatrixAdaptor.h"
#include "AMP/AMP_TPLs.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/data/hypre/HypreCSRPolicy.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/memory.h"

#include <numeric>

extern "C" {
#include "HYPRE_utilities.h"
#include "_hypre_IJ_mv.h"
#include "_hypre_parcsr_mv.h"
}

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
        csrData->getOffDiagColumnMap( d_colMap );

        AMP_INSIST( nnz_d && cols_d && cols_loc_d && coeffs_d,
                    "diagonal block layout cannot be NULL" );

        initializeHypreMatrix( csrData->getMemoryLocation(),
                               firstRow,
                               lastRow,
                               csrData->numberOfNonZerosDiag(),
                               nnz_d,
                               cols_d,
                               cols_loc_d,
                               coeffs_d,
                               csrData->numberOfNonZerosOffDiag(),
                               nnz_od,
                               cols_od,
                               cols_loc_od,
                               coeffs_od,
                               d_colMap.size(),
                               d_colMap.data() );

    } else {

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

HypreMatrixAdaptor::~HypreMatrixAdaptor()
{
    hypre_ParCSRMatrix *par_matrix = static_cast<hypre_ParCSRMatrix *>( d_matrix->object );
    par_matrix->col_map_offd       = nullptr;
    // Now the standard IJMatrixDestroy can be called
    HYPRE_IJMatrixDestroy( d_matrix );
}

void HypreMatrixAdaptor::initializeHypreMatrix( AMP::Utilities::MemoryType mem_loc,
                                                HYPRE_BigInt first_row,
                                                HYPRE_BigInt last_row,
                                                HYPRE_BigInt nnz_total_d,
                                                HYPRE_Int *nnz_per_row_d,
                                                HYPRE_BigInt *csr_bja_d,
                                                HYPRE_Int *csr_lja_d,
                                                HYPRE_Real *csr_aa_d,
                                                HYPRE_BigInt nnz_total_od,
                                                HYPRE_Int *nnz_per_row_od,
                                                HYPRE_BigInt *csr_bja_od,
                                                HYPRE_Int *csr_lja_od,
                                                HYPRE_Real *csr_aa_od,
                                                HYPRE_BigInt csr_col_map_size,
                                                HYPRE_BigInt *csr_col_map_offd )
{
    if ( mem_loc == AMP::Utilities::MemoryType::host ) {
        HYPRE_SetMemoryLocation( HYPRE_MEMORY_HOST );
    } else {
        AMP_ERROR( "Non-host memory not yet supported in HypreMatrixAdaptor" );
    }

    const auto nrows = last_row - first_row + 1;

    // Manually create ParCSR and fill fields as needed
    // Roughly based on hypre_IJMatrixInitializeParCSR_v2 from IJMatrix_parcsr.c
    //   and the various functions that it calls
    hypre_IJMatrixCreateParCSR( d_matrix );
    hypre_ParCSRMatrix *par_matrix = static_cast<hypre_ParCSRMatrix *>( d_matrix->object );
    hypre_CSRMatrix *diag          = par_matrix->diag;
    hypre_CSRMatrix *off_diag      = par_matrix->offd;

    // Filling the contents manually should remove any need for aux matrix
    hypre_AuxParCSRMatrix *aux_mat = static_cast<hypre_AuxParCSRMatrix *>( d_matrix->translator );
    aux_mat->need_aux              = 0;

    // Verify that Hypre CSRMatrices are on host memory
    AMP_INSIST( diag->memory_location == HYPRE_MEMORY_HOST &&
                    off_diag->memory_location == HYPRE_MEMORY_HOST,
                "Hypre matrices need to be on host memory for adaptor to work" );

    // Verify that diag and off_diag are "empty"
    AMP_INSIST( diag->num_nonzeros == 0 && off_diag->num_nonzeros == 0,
                "Hypre (off)diag matrix has nonzeros but shouldn't" );

    // Hypre always frees the hypre_CSRMatrix->i and hypre_CSRMatrix->rownnz
    // fields regardless of ->owns_data. Calling matrix initialize will let
    // hypre do those allocations. ->big_j, ->j, and ->data should not get
    // allocated since ->num_nonzeros == 0
    hypre_CSRMatrixInitialize( diag );
    hypre_CSRMatrixInitialize( off_diag );

    // Fill in the ->i and ->rownnz fields of diag and off_diag
    diag->i[0]     = 0;
    off_diag->i[0] = 0;
    for ( HYPRE_BigInt n = 0; n < nrows; ++n ) {
        diag->i[n + 1]     = diag->i[n] + nnz_per_row_d[n];
        off_diag->i[n + 1] = off_diag->i[n] + nnz_per_row_od[n];
    }

    // This is where we tell hypre to stop owning any data
    hypre_CSRMatrixSetDataOwner( diag, 0 );
    hypre_CSRMatrixSetDataOwner( off_diag, 0 );

    // Now set diag/off_diag members to point at our data
    diag->big_j = csr_bja_d;
    diag->data  = csr_aa_d;
    diag->j     = csr_lja_d;

    off_diag->big_j = csr_bja_od;
    off_diag->data  = csr_aa_od;
    off_diag->j     = csr_lja_od;

    // Update metadata fields to match what we've inserted
    diag->num_nonzeros     = nnz_total_d;
    off_diag->num_nonzeros = nnz_total_od;

    // permute diag->j and ->data so that the diagonal element is first
    // IJMatrix_parcsr.c in function
    for ( HYPRE_Int i = 0; i < static_cast<HYPRE_Int>( nrows ); ++i ) {
        auto j0 = diag->i[i];
        for ( HYPRE_Int j = j0; j < diag->i[i + 1]; ++j ) {
            // Hypre does not permute diag->big_j but we need to for consistency
            if ( diag->j[j] == i ) {
                auto dTemp      = diag->data[j0];
                auto bjTmp      = diag->big_j[j0];
                diag->data[j0]  = diag->data[j];
                diag->data[j]   = dTemp;
                diag->big_j[j0] = diag->big_j[j];
                diag->big_j[j]  = bjTmp;
                diag->j[j]      = diag->j[j0];
                diag->j[j0]     = i;
                break;
            }
        }
    }

    // Set colmap inside ParCSR and flag that assembly is already done
    // See destructor above regarding ownership of this field
    par_matrix->col_map_offd = csr_col_map_offd;
    off_diag->num_cols       = csr_col_map_size;

    // Update ->rownnz fields, note that we don't own these
    hypre_CSRMatrixSetRownnz( diag );
    hypre_CSRMatrixSetRownnz( off_diag );

    // set assemble flag to indicate that we are done
    d_matrix->assemble_flag = 1;
}

} // namespace AMP::LinearAlgebra
