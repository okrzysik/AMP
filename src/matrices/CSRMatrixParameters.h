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

    // The diagonal and off-diagonal blocks need all the same parameters
    // Like in CSRMatrixData use a nested class to pack all this away
    struct CSRLocalMatrixParameters {
        // No bare constructor, only initializer lists and default copy/moves
        CSRLocalMatrixParameters() = delete;

        lidx_t *d_row_starts;
        gidx_t *d_cols;
        scalar_t *d_coeffs;
    };

    CSRMatrixParameters() = delete;

    /** \brief Constructor
     * \param[in] first_row     Index for first row
     * \param[in] last_row      Index for last row
     * \param[in] first_col     Index for first col
     * \param[in] last_col      Index for last col
     * \param[in] diag          Parameters for diag block
     * \param[in] off_diag      Parameters for offd block
     * \param[in] comm          Communicator for the matrix
     */
    explicit CSRMatrixParameters( gidx_t first_row,
                                  gidx_t last_row,
                                  gidx_t first_col,
                                  gidx_t last_col,
                                  const CSRLocalMatrixParameters &diag,
                                  const CSRLocalMatrixParameters &off_diag,
                                  const AMP_MPI &comm )
        : MatrixParametersBase( comm ),
          d_first_row( first_row ),
          d_last_row( last_row ),
          d_first_col( first_col ),
          d_last_col( last_col ),
          d_diag( diag ),
          d_off_diag( off_diag )
    {
    }

    //! Destructor
    virtual ~CSRMatrixParameters() = default;

    // Bulk information
    gidx_t d_first_row;
    gidx_t d_last_row;
    gidx_t d_first_col;
    gidx_t d_last_col;
    // Blockwise information
    CSRLocalMatrixParameters d_diag, d_off_diag;
};
} // namespace AMP::LinearAlgebra

#endif
