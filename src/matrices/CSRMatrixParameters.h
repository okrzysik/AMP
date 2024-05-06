#ifndef included_AMP_CSRMatrixParameters
#define included_AMP_CSRMatrixParameters

#include "AMP/matrices/MatrixParametersBase.h"

namespace AMP::LinearAlgebra {


/** \class CSRMatrixParameters
 * \brief  A class used to hold basic parameters for a matrix
 */
template<typename CSRPolicy>
class CSRMatrixParameters : public MatrixParametersBase
{
public:
    using gidx_t   = typename CSRPolicy::gidx_t;
    using lidx_t   = typename CSRPolicy::lidx_t;
    using scalar_t = typename CSRPolicy::scalar_t;

    CSRMatrixParameters() = delete;

    /** \brief Constructor
     * \param[in] comm     Communicator for the matrix
     */
    explicit CSRMatrixParameters( gidx_t first_row,
                                  gidx_t last_row,
                                  lidx_t *nnz_per_row_diag,
                                  gidx_t *cols_diag,
                                  lidx_t *cols_loc_diag,
                                  scalar_t *coeffs_diag,
                                  lidx_t *nnz_per_row_odiag,
                                  gidx_t *cols_odiag,
                                  lidx_t *cols_loc_odiag,
                                  scalar_t *coeffs_odiag,
                                  const AMP_MPI &comm )
        : MatrixParametersBase( comm ),
          d_is_square( true ),
          d_first_row( first_row ),
          d_last_row( last_row ),
          d_first_col( first_row ),
          d_last_col( last_row ),
          d_nnz_per_row_diag( nnz_per_row_diag ),
          d_cols_diag( cols_diag ),
          d_cols_loc_diag( cols_loc_diag ),
          d_coeffs_diag( coeffs_diag ),
          d_nnz_per_row_odiag( nnz_per_row_odiag ),
          d_cols_odiag( cols_odiag ),
          d_cols_loc_odiag( cols_loc_odiag ),
          d_coeffs_odiag( coeffs_odiag )
    {
    }

    explicit CSRMatrixParameters( bool is_square,
                                  gidx_t first_row,
                                  gidx_t last_row,
                                  gidx_t first_col,
                                  gidx_t last_col,
                                  lidx_t *nnz_per_row_diag,
                                  gidx_t *cols_diag,
                                  lidx_t *cols_loc_diag,
                                  scalar_t *coeffs_diag,
                                  lidx_t *nnz_per_row_odiag,
                                  gidx_t *cols_odiag,
                                  lidx_t *cols_loc_odiag,
                                  scalar_t *coeffs_odiag,
                                  const AMP_MPI &comm )
        : MatrixParametersBase( comm ),
          d_is_square( is_square ),
          d_first_row( first_row ),
          d_last_row( last_row ),
          d_first_col( first_col ),
          d_last_col( last_col ),
          d_nnz_per_row_diag( nnz_per_row_diag ),
          d_cols_diag( cols_diag ),
          d_cols_loc_diag( cols_loc_diag ),
          d_coeffs_diag( coeffs_diag ),
          d_nnz_per_row_odiag( nnz_per_row_odiag ),
          d_cols_odiag( cols_odiag ),
          d_cols_loc_odiag( cols_loc_odiag ),
          d_coeffs_odiag( coeffs_odiag )
    {
    }

    //! Deconstructor
    virtual ~CSRMatrixParameters() = default;
    // Overall information
    bool d_is_square;
    gidx_t d_first_row;
    gidx_t d_last_row;
    gidx_t d_first_col;
    gidx_t d_last_col;
    // Diagonal block info
    lidx_t *d_nnz_per_row_diag;
    gidx_t *d_cols_diag;
    lidx_t *d_cols_loc_diag;
    scalar_t *d_coeffs_diag;
    // Off-Diagonal block info
    lidx_t *d_nnz_per_row_odiag;
    gidx_t *d_cols_odiag;
    lidx_t *d_cols_loc_odiag;
    scalar_t *d_coeffs_odiag;
};
} // namespace AMP::LinearAlgebra

#endif
