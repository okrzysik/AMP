#ifndef included_AMP_CSRSerialMatrixParameters
#define included_AMP_CSRSerialMatrixParameters

#include "AMP/matrices/MatrixParametersBase.h"

namespace AMP::LinearAlgebra {


/** \class CSRSerialMatrixParameters
 * \brief  A class used to hold basic parameters for a matrix
 */
template<typename CSRPolicy>
class CSRSerialMatrixParameters : public MatrixParametersBase
{
public:
    using gidx_t   = typename CSRPolicy::gidx_t;
    using lidx_t   = typename CSRPolicy::lidx_t;
    using scalar_t = typename CSRPolicy::scalar_t;

    CSRSerialMatrixParameters() = delete;

    explicit CSRSerialMatrixParameters( bool is_diag,
					bool is_empty,
					gidx_t first_row,
					gidx_t last_row,
					gidx_t first_col,
					gidx_t last_col,
					lidx_t *nnz_per_row,
					gidx_t *cols,
					lidx_t *cols_loc,
					scalar_t *coeffs,
					const AMP_MPI &comm )
        : MatrixParametersBase( comm ),
	  d_is_diag( is_diag ),
	  d_is_empty( is_empty ),
	  d_first_row( first_row ),
          d_last_row( last_row ),
          d_first_col( first_col ),
          d_last_col( last_col ),
          d_nnz_per_row( nnz_per_row ),
          d_cols( cols ),
          d_cols_loc( cols_loc ),
          d_coeffs( coeffs )
    {
    }

    //! Deconstructor
    virtual ~CSRSerialMatrixParameters() = default;

    bool d_is_diag;
    bool d_is_empty;
    gidx_t d_first_row;
    gidx_t d_last_row;
    gidx_t d_first_col;
    gidx_t d_last_col;
    lidx_t *d_nnz_per_row;
    gidx_t *d_cols;
    lidx_t *d_cols_loc;
    scalar_t *d_coeffs;
};
} // namespace AMP::LinearAlgebra

#endif
