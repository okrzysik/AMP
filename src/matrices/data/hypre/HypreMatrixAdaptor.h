#ifndef INCLUDED_HYPRE_MATRIX_ADAPTOR_H
#define INCLUDED_HYPRE_MATRIX_ADAPTOR_H

#include "AMP/matrices/data/MatrixData.h"

extern "C" {
#include "_hypre_IJ_mv.h"
#include "_hypre_parcsr_mv.h"

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_utilities.h"
}

#include <vector>

namespace AMP::LinearAlgebra {

/** \class HypreMatrixAdaptor
 * \brief  Wrapper class for constructing a hypre IJ matrix for use with solvers
 * \details It is initialized based on CSR data.
 */
class HypreMatrixAdaptor
{

public:
    HypreMatrixAdaptor() = delete;
    ~HypreMatrixAdaptor();

    //!! Main constructor using a MatrixData object
    HypreMatrixAdaptor( std::shared_ptr<MatrixData> );

    //! Returns the handle to the HYPRE IJMatrix
    HYPRE_IJMatrix getHypreMatrix( void ) { return d_matrix; }

private:
    //! Main internal routine for initializing the matrix
    void initializeHypreMatrix( AMP::Utilities::MemoryType,
                                HYPRE_BigInt first_row,
                                HYPRE_BigInt last_row,
                                bool has_off_diag,
                                HYPRE_BigInt nnz_total_d,
                                HYPRE_Int *csr_ia_d,
                                HYPRE_BigInt *csr_bja_d,
                                HYPRE_Int *csr_lja_d,
                                HYPRE_Real *csr_aa_d,
                                HYPRE_BigInt nnz_total_od,
                                HYPRE_Int *csr_ia_od,
                                HYPRE_BigInt *csr_bja_od,
                                HYPRE_Int *csr_lja_od,
                                HYPRE_Real *csr_aa_od,
                                HYPRE_BigInt csr_col_map_size,
                                HYPRE_BigInt *csr_col_map );

    //! hypre IJ matrix that this class wraps
    HYPRE_IJMatrix d_matrix;

    std::vector<HYPRE_BigInt> d_colMap;
};

} // namespace AMP::LinearAlgebra

#endif // INCLUDED_HYPRE_MATRIX_ADAPTOR_H
