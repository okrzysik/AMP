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


template void HypreMatrixAdaptor::initializeHypreMatrix<
    CSRMatrixData<HypreCSRPolicy, AMP::HostAllocator<int>>>(
    std::shared_ptr<CSRMatrixData<HypreCSRPolicy, AMP::HostAllocator<int>>> );

#ifdef USE_DEVICE
template void HypreMatrixAdaptor::initializeHypreMatrix<
    CSRMatrixData<HypreCSRPolicy, AMP::ManagedAllocator<int>>>(
    std::shared_ptr<CSRMatrixData<HypreCSRPolicy, AMP::ManagedAllocator<int>>> );

template void HypreMatrixAdaptor::initializeHypreMatrix<
    CSRMatrixData<HypreCSRPolicy, AMP::DeviceAllocator<int>>>(
    std::shared_ptr<CSRMatrixData<HypreCSRPolicy, AMP::DeviceAllocator<int>>> );
#endif

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

    // Attempt dynamic pointer casts to supported types
    // Policy must match HypreCSRPolicy
    // need to match supported allocators depending on device support
    auto csrDataHost =
        std::dynamic_pointer_cast<CSRMatrixData<HypreCSRPolicy, AMP::HostAllocator<int>>>(
            matrixData );

#ifdef USE_DEVICE
    auto csrDataManaged =
        std::dynamic_pointer_cast<CSRMatrixData<HypreCSRPolicy, AMP::ManagedAllocator<int>>>(
            matrixData );
    auto csrDataDevice =
        std::dynamic_pointer_cast<CSRMatrixData<HypreCSRPolicy, AMP::DeviceAllocator<int>>>(
            matrixData );
#else
    // Just default out these to nullptrs to make logic below simpler
    decltype( csrDataHost ) csrDataManaged = nullptr;
    decltype( csrDataHost ) csrDataDevice  = nullptr;
#endif

    if ( csrDataHost ) {
        initializeHypreMatrix( csrDataHost );
    } else if ( csrDataManaged && false ) {
        initializeHypreMatrix( csrDataManaged );
    } else if ( csrDataDevice && false ) {
        initializeHypreMatrix( csrDataDevice );
    } else {

        HYPRE_SetMemoryLocation( HYPRE_MEMORY_HOST );
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

template<class csr_data_type>
void HypreMatrixAdaptor::initializeHypreMatrix( std::shared_ptr<csr_data_type> csrData )
{
    // extract fields from csrData
    HYPRE_BigInt first_row = static_cast<HYPRE_BigInt>( csrData->beginRow() );
    HYPRE_BigInt last_row  = static_cast<HYPRE_BigInt>( csrData->endRow() - 1 );
    csrData->getOffDiagColumnMap( d_colMap );
    HYPRE_BigInt nnz_total_d  = static_cast<HYPRE_BigInt>( csrData->numberOfNonZerosDiag() );
    HYPRE_BigInt nnz_total_od = static_cast<HYPRE_BigInt>( csrData->numberOfNonZerosOffDiag() );
    auto [nnz_d, cols_d, cols_loc_d, coeffs_d]     = csrData->getCSRDiagData();
    auto [nnz_od, cols_od, cols_loc_od, coeffs_od] = csrData->getCSROffDiagData();

    AMP_INSIST( nnz_d && cols_d && cols_loc_d && coeffs_d, "diagonal block layout cannot be NULL" );

    if ( csrData->getMemoryLocation() == AMP::Utilities::MemoryType::host ) {
        HYPRE_SetMemoryLocation( HYPRE_MEMORY_HOST );
    } else if ( csrData->getMemoryLocation() > AMP::Utilities::MemoryType::host ) {
#ifdef USE_DEVICE
        HYPRE_SetMemoryLocation( HYPRE_MEMORY_DEVICE );
        AMP_ERROR( "Non-host memory not yet supported in HypreMatrixAdaptor" );
#else
        AMP_ERROR( "Non-host memory not yet supported in HypreMatrixAdaptor" );
#endif
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
        diag->i[n + 1]     = diag->i[n] + nnz_d[n];
        off_diag->i[n + 1] = off_diag->i[n] + nnz_od[n];
    }

    // This is where we tell hypre to stop owning any data
    hypre_CSRMatrixSetDataOwner( diag, 0 );
    hypre_CSRMatrixSetDataOwner( off_diag, 0 );

    // Now set diag/off_diag members to point at our data
    diag->big_j = cols_d;
    diag->data  = coeffs_d;
    diag->j     = cols_loc_d;

    off_diag->big_j = cols_od;
    off_diag->data  = coeffs_od;
    off_diag->j     = cols_loc_od;

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
    par_matrix->col_map_offd = d_colMap.data();
    off_diag->num_cols       = d_colMap.size();

    // Update ->rownnz fields, note that we don't own these
    hypre_CSRMatrixSetRownnz( diag );
    hypre_CSRMatrixSetRownnz( off_diag );

    // set assemble flag to indicate that we are done
    d_matrix->assemble_flag = 1;
}

} // namespace AMP::LinearAlgebra
