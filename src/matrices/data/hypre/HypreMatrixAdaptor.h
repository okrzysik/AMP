#ifndef INCLUDED_HYPRE_MATRIX_ADAPTOR_H
#define INCLUDED_HYPRE_MATRIX_ADAPTOR_H

#include "AMP/matrices/data/MatrixData.h"

extern "C" {
#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_utilities.h"
}

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
    void initializeHypreMatrix( HYPRE_BigInt first_row,
                                HYPRE_BigInt last_row,
                                bool has_off_diag,
                                HYPRE_Int *const csr_ia_d,
                                HYPRE_BigInt const *const csr_ja_d,
                                HYPRE_Real const *const csr_aa_d,
                                HYPRE_Int *const csr_ia_od,
                                HYPRE_BigInt const *const csr_ja_od,
                                HYPRE_Real const *const csr_aa_od );

    //! hypre IJ matrix that this class wraps
    HYPRE_IJMatrix d_matrix;
};

} // namespace AMP::LinearAlgebra

#endif // INCLUDED_HYPRE_MATRIX_ADAPTOR_H
