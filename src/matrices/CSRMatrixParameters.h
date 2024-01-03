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
                                  lidx_t *nnz_per_row,
                                  gidx_t const *const cols,
                                  scalar_t const *const coeffs,
                                  const AMP_MPI &comm )
        : MatrixParametersBase( comm ),
          d_is_square( true ),
          d_first_row( first_row ),
          d_last_row( last_row ),
          d_first_col( first_row ),
          d_last_col( last_row ),
          d_nnz_per_row( nnz_per_row ),
          d_cols( cols ),
          d_coeffs( coeffs )
    {
    }

    explicit CSRMatrixParameters( bool is_square,
                                  gidx_t first_row,
                                  gidx_t last_row,
                                  gidx_t first_col,
                                  gidx_t last_col,
                                  lidx_t *nnz_per_row,
                                  gidx_t const *const cols,
                                  scalar_t const *const coeffs,
                                  const AMP_MPI &comm )
        : MatrixParametersBase( comm ),
          d_is_square( is_square ),
          d_first_row( first_row ),
          d_last_row( last_row ),
          d_first_col( first_col ),
          d_last_col( last_col ),
          d_nnz_per_row( nnz_per_row ),
          d_cols( cols ),
          d_coeffs( coeffs )
    {
    }

    //! Deconstructor
    virtual ~CSRMatrixParameters() = default;

    bool d_is_square;
    gidx_t d_first_row;
    gidx_t d_last_row;
    gidx_t d_first_col;
    gidx_t d_last_col;
    lidx_t *d_nnz_per_row;
    gidx_t const *d_cols;
    scalar_t const *d_coeffs;
};
} // namespace AMP::LinearAlgebra

#endif
