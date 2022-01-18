#ifndef included_AMP_ManagedMatrix
#define included_AMP_ManagedMatrix

#include "Matrix.h"


namespace AMP::LinearAlgebra {


/** \class ManagedMatrix
 * \brief  A class that indicates AMP manages the data is some fashion
 * \details  This class is similar to the ManagedVector class in that
 * it indicates AMP manages the memory in some fashion.  ManagedMatrix
 * will provide an interface sufficient to implement a CRS matrix.
 */
class ManagedMatrix : virtual public Matrix
{
protected:
    ManagedMatrix();

public:
    explicit ManagedMatrix( std::shared_ptr<ManagedMatrixParameters> p );

    /** \brief  Create values in a matrix
     * \param[in]  row       the row ids
     * \param[in]  cols      the column ids
     * \details  A managed matrix may require a separate step of
     * creating non-zeros to be packed later.  Rather than using
     * set, this method is used.  Once all createValuesByGlobalID
     * are finished, call fillComplete().
     */
    virtual void createValuesByGlobalID( size_t row, const std::vector<size_t> &cols ) = 0;

    /** \brief  All createValuesByGlobalID have completed.
     */
    virtual void fillComplete() = 0;

protected:
    //!  Parameters used to construct the matrix
    std::shared_ptr<ManagedMatrixParameters> d_pParameters;
};

} // namespace AMP::LinearAlgebra

#include "ManagedMatrix.inline.h"
#endif
