#ifndef included_AMP_CSRMatrix
#define included_AMP_CSRMatrix

#include "AMP/matrices/Matrix.h"
#include "AMP/vectors/Vector.h"
#include <functional>
#include <memory>

namespace AMP::LinearAlgebra {

/** \class CSRMatrix
 * \brief  An concrete class for dealing with dense serial matrices
 * \details  This is a concrete class that stores a dense local matrix.
 *    This is not a distributed matrix and requires that the comm is AMP_COMM_SELF.
 */
template<typename Policy, typename Allocator = std::allocator<int>>
class CSRMatrix : public Matrix
{
public:
    CSRMatrix() = delete;

    /** \brief Constructor
     * \param[in] params  Description of the matrix
     */
    explicit CSRMatrix( std::shared_ptr<MatrixParametersBase> params );

    /** \brief Constructor
     * \param[in] params  MatrixData object associated with matrix
     */
    explicit CSRMatrix( std::shared_ptr<MatrixData> data );

    CSRMatrix( const CSRMatrix & ) = delete;

    CSRMatrix &operator=( const CSRMatrix & ) = delete;

    /** \brief Destructor
     */
    virtual ~CSRMatrix();

    //! Return the type of the matrix
    virtual std::string type() const override { return "CSRMatrix"; }

    /** \brief  Return a new matrix that is the transpose of this one
     * \return  A copy of this matrix transposed.
     */
    shared_ptr transpose() const override;

    /** \brief  Return a matrix with the same non-zero and distributed structure
     * \return  The new matrix
     */
    shared_ptr clone() const override;

    /** \brief  Extract the diagonal from a matrix
     * \param[in]  buf  An optional vector to use as a buffer
     * \return  A vector of the diagonal values
     */
    Vector::shared_ptr
    extractDiagonal( Vector::shared_ptr buf = Vector::shared_ptr() ) const override;

    /** \brief Get a right vector( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$\mathbf{x}\f$ is a right
     * vector )
     * \return  A newly created right vector
     */
    Vector::shared_ptr getRightVector() const override;

    /** \brief Get a left vector( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$\mathbf{y}\f$ is a left
     * vector )
     * \return  A newly created left vector
     */
    Vector::shared_ptr getLeftVector() const override;

protected:
    /** \brief  Multiply two matrices and store in a third
     * \param[in]  other_op  The other matrix to multiply
     * \param[out] result  The matrix to store the result
     */
    void multiply( shared_ptr other_op, shared_ptr &result ) override;
};
} // namespace AMP::LinearAlgebra


#endif
