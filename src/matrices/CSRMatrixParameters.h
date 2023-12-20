#ifndef included_AMP_CSRMatrixParameters
#define included_AMP_CSRMatrixParameters

#include "AMP/matrices/MatrixParametersBase.h"

namespace AMP::LinearAlgebra {


/** \class CSRMatrixParameters
 * \brief  A class used to hold basic parameters for a matrix
 */
class CSRMatrixParameters : public MatrixParametersBase
{
public:
    CSRMatrixParameters() = delete;

    /** \brief Constructor
     * \param[in] comm     Communicator for the matrix
     */
    explicit CSRMatrixParameters( size_t first_row,
                                  size_t last_row,
                                  size_t const *const nnz_per_row,
                                  size_t const *const cols,
                                  double const *const coeffs,
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
                                  size_t first_row,
                                  size_t last_row,
                                  size_t first_col,
                                  size_t last_col,
                                  size_t const *const nnz_per_row,
                                  size_t const *const cols,
                                  double const *const coeffs,
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
    size_t d_first_row;
    size_t d_last_row;
    size_t d_first_col;
    size_t d_last_col;
    size_t const *d_nnz_per_row;
    size_t const *d_cols;
    double const *d_coeffs;
};
} // namespace AMP::LinearAlgebra

#endif
