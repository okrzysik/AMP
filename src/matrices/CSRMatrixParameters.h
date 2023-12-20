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
                                  int32_t const *const nnz_per_row,
                                  int64_t const *const cols,
                                  double const *const coeffs,
                                  const AMP_MPI &comm )
        : MatrixParametersBase( comm ),
          d_first_row( first_row ),
          d_last_row( last_row ),
          d_nnz_per_row( nnz_per_row ),
          d_cols( cols ),
          d_coeffs( coeffs )
    {
    }

    //! Deconstructor
    virtual ~CSRMatrixParameters() = default;

    size_t d_first_row;
    size_t d_last_row;
    int32_t const *d_nnz_per_row;
    int64_t const *d_cols;
    double const *d_coeffs;
};
} // namespace AMP::LinearAlgebra

#endif
